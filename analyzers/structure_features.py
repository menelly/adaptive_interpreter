#!/usr/bin/env python3
"""
Structure features — numpy-native AlphaFold per-residue geometry.

Mechanism-first variant analysis needs to know the *physical context* of a
residue: is it buried in the fold or on the surface? in a rigid helix/sheet or
a flexible loop? who does it touch in 3D (salt-bridge latches, interfaces)? how
confident/ordered is the local structure (pLDDT)?

We extract all of that straight from the AlphaFold human-proteome PDBs with
nothing but numpy — deliberately NO Biopython dependency (it isn't installed in
the shared codex venv, and the old structural path silently fell back to
defaults because of that). Pure ATOM-record parsing + geometry.

No hardcoded genes. No frequency. No conservation. Just biophysics from the
single-chain AlphaFold model.

Cache lives on the big D mount (arcana + the main drive are tight):
    /mnt/win-d/Ace/adaptiveinterpreter_cache/structure_features/{uniprot}.npz
"""
from __future__ import annotations

import os
import gzip
import logging
from typing import Dict, List, Optional

import numpy as np

logger = logging.getLogger(__name__)

AF_STRUCT_DIR = "/mnt/arcana/alphafold_human/structures"
CACHE_DIR = "/mnt/win-d/Ace/adaptiveinterpreter_cache/structure_features"

# Burial: number of CB neighbours within this radius (coordination number).
BURIAL_RADIUS = 10.0
# Contact map: CB-CB closer than this are "in contact".
CONTACT_RADIUS = 8.0
# Salt bridge: oppositely-charged side-chain centroids within this distance.
SALT_BRIDGE_RADIUS = 5.0
# Coordination-number thresholds (calibrated against typical globular CN ~ 8-22).
SURFACE_CN = 11.0   # at or below -> surface-exposed
BURIED_CN = 18.0    # at or above -> buried core

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

POSITIVE = set("KR")          # His left out: weakly/conditionally charged
NEGATIVE = set("DE")


def af_path(uniprot_id: str) -> Optional[str]:
    """Return the AlphaFold PDB path for a UniProt accession (fragment 1)."""
    p = os.path.join(AF_STRUCT_DIR, f"AF-{uniprot_id}-F1-model_v4.pdb.gz")
    return p if os.path.exists(p) else None


def _dihedral(p0, p1, p2, p3) -> float:
    """Signed dihedral angle (degrees) defined by four points."""
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / (np.linalg.norm(b1) + 1e-9)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    x = np.dot(v, w)
    y = np.dot(np.cross(b1n, v), w)
    return float(np.degrees(np.arctan2(y, x)))


class StructureFeatures:
    """Per-residue structural features for one protein, lazily cached.

    Residues are addressed by 1-based author sequence position (resSeq), which
    matches HGVS protein numbering for the canonical AlphaFold isoform.
    """

    def __init__(self, uniprot_id: str):
        self.uniprot_id = uniprot_id
        self.ok = False
        self.positions: np.ndarray = np.array([], dtype=int)
        self.aa: Dict[int, str] = {}
        self.plddt: Dict[int, float] = {}
        self.cn: Dict[int, float] = {}              # coordination number (burial)
        self.ss: Dict[int, str] = {}                # 'H' helix / 'E' sheet / 'C' coil
        self._contacts: Dict[int, List[int]] = {}
        self._sc_centroid: Dict[int, np.ndarray] = {}
        self._load()

    # ---- public per-residue accessors -------------------------------------

    def has(self, pos: int) -> bool:
        return self.ok and pos in self.aa

    def burial(self, pos: int) -> str:
        """'surface' | 'buried' | 'intermediate' for a residue."""
        cn = self.cn.get(pos)
        if cn is None:
            return "unknown"
        if cn <= SURFACE_CN:
            return "surface"
        if cn >= BURIED_CN:
            return "buried"
        return "intermediate"

    def is_disordered(self, pos: int) -> bool:
        """Low pLDDT => flexible/disordered (where degrons & regulatory sites live)."""
        return self.plddt.get(pos, 100.0) < 70.0

    def contacts(self, pos: int) -> List[int]:
        return self._contacts.get(pos, [])

    def salt_bridge_partners(self, pos: int) -> List[int]:
        """Residues whose side chain forms a salt bridge with `pos` (opposite charge, <5A)."""
        a = self.aa.get(pos)
        if a not in POSITIVE and a not in NEGATIVE:
            return []
        want = NEGATIVE if a in POSITIVE else POSITIVE
        c0 = self._sc_centroid.get(pos)
        if c0 is None:
            return []
        out = []
        for j in self._contacts.get(pos, []):
            if self.aa.get(j) in want:
                cj = self._sc_centroid.get(j)
                if cj is not None and np.linalg.norm(c0 - cj) <= SALT_BRIDGE_RADIUS:
                    out.append(j)
        return out

    def summary(self, pos: int) -> Dict:
        return {
            "uniprot": self.uniprot_id,
            "aa": self.aa.get(pos),
            "plddt": round(self.plddt.get(pos, float("nan")), 1),
            "coordination_number": round(self.cn.get(pos, float("nan")), 1),
            "burial": self.burial(pos),
            "secondary_structure": self.ss.get(pos),
            "disordered": self.is_disordered(pos),
            "n_contacts": len(self._contacts.get(pos, [])),
            "salt_bridge_partners": self.salt_bridge_partners(pos),
        }

    # ---- loading / caching ------------------------------------------------

    def _cache_file(self) -> str:
        return os.path.join(CACHE_DIR, f"{self.uniprot_id}.npz")

    def _load(self):
        if self._load_cache():
            self.ok = True
            return
        path = af_path(self.uniprot_id)
        if not path:
            logger.warning(f"No AlphaFold structure for {self.uniprot_id}")
            return
        try:
            self._parse_and_compute(path)
            self.ok = True
            self._save_cache()
        except Exception as e:  # never let structure failure break the cascade
            logger.warning(f"Structure parse failed for {self.uniprot_id}: {e}")
            self.ok = False

    def _parse_and_compute(self, path: str):
        # backbone atoms for dihedrals + side-chain atoms for centroids
        bb: Dict[int, Dict[str, np.ndarray]] = {}   # pos -> {'N','CA','C'}
        cb: Dict[int, np.ndarray] = {}              # pos -> CB (CA for Gly)
        sc_atoms: Dict[int, List[np.ndarray]] = {}  # pos -> side-chain coords
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if not line.startswith("ATOM"):
                    continue
                atom = line[12:16].strip()
                alt = line[16]
                if alt not in (" ", "A"):
                    continue
                resn = line[17:20].strip()
                aa1 = THREE_TO_ONE.get(resn)
                if aa1 is None:
                    continue
                pos = int(line[22:26])
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                self.aa.setdefault(pos, aa1)
                if atom in ("N", "CA", "C"):
                    bb.setdefault(pos, {})[atom] = xyz
                    if atom == "CA":
                        self.plddt[pos] = float(line[60:66])
                if atom == "CB" or (atom == "CA" and aa1 == "G"):
                    cb[pos] = xyz
                if atom not in ("N", "CA", "C", "O"):
                    sc_atoms.setdefault(pos, []).append(xyz)

        self.positions = np.array(sorted(self.aa.keys()), dtype=int)
        # side-chain centroid (fallback to CB/CA) for salt-bridge geometry
        for pos in self.positions:
            if sc_atoms.get(pos):
                self._sc_centroid[pos] = np.mean(sc_atoms[pos], axis=0)
            elif pos in cb:
                self._sc_centroid[pos] = cb[pos]

        # coordination number + contact map from CB coords (vectorised)
        cb_pos = [p for p in self.positions if p in cb]
        coords = np.array([cb[p] for p in cb_pos])
        if len(coords):
            d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
            np.fill_diagonal(d, np.inf)
            for i, p in enumerate(cb_pos):
                self.cn[p] = float(np.sum(d[i] <= BURIAL_RADIUS))
                self._contacts[p] = [cb_pos[j] for j in np.where(d[i] <= CONTACT_RADIUS)[0]]

        # secondary structure from backbone phi/psi
        self._assign_secondary_structure(bb)

    def _assign_secondary_structure(self, bb: Dict[int, Dict[str, np.ndarray]]):
        for pos in self.positions:
            ss = "C"
            prev, cur, nxt = bb.get(pos - 1), bb.get(pos), bb.get(pos + 1)
            if cur and prev and nxt and all(k in cur for k in ("N", "CA", "C")) \
                    and "C" in prev and "N" in nxt:
                phi = _dihedral(prev["C"], cur["N"], cur["CA"], cur["C"])
                psi = _dihedral(cur["N"], cur["CA"], cur["C"], nxt["N"])
                if -150 <= phi <= -30 and -90 <= psi <= 30:
                    ss = "H"   # alpha-helix basin
                elif -180 <= phi <= -90 and (90 <= psi <= 180 or -180 <= psi <= -150):
                    ss = "E"   # beta-strand basin
            self.ss[pos] = ss

    def _save_cache(self):
        try:
            os.makedirs(CACHE_DIR, exist_ok=True)
            np.savez_compressed(
                self._cache_file(),
                positions=self.positions,
                aa=np.array([self.aa[p] for p in self.positions]),
                plddt=np.array([self.plddt.get(p, np.nan) for p in self.positions]),
                cn=np.array([self.cn.get(p, np.nan) for p in self.positions]),
                ss=np.array([self.ss.get(p, "C") for p in self.positions]),
                sc_cx=np.array([self._sc_centroid.get(p, [np.nan] * 3) for p in self.positions]),
                contacts=np.array([",".join(map(str, self._contacts.get(p, []))) for p in self.positions]),
            )
        except Exception as e:
            logger.warning(f"Could not cache {self.uniprot_id}: {e}")

    def _load_cache(self) -> bool:
        f = self._cache_file()
        if not os.path.exists(f):
            return False
        try:
            z = np.load(f, allow_pickle=False)
            self.positions = z["positions"].astype(int)
            aa, pl, cn, ss = z["aa"], z["plddt"], z["cn"], z["ss"]
            sc, ct = z["sc_cx"], z["contacts"]
            for i, p in enumerate(self.positions):
                p = int(p)
                self.aa[p] = str(aa[i])
                self.plddt[p] = float(pl[i])
                self.cn[p] = float(cn[i])
                self.ss[p] = str(ss[i])
                self._sc_centroid[p] = np.array(sc[i], dtype=float)
                self._contacts[p] = [int(x) for x in str(ct[i]).split(",") if x]
            return True
        except Exception as e:
            logger.warning(f"Cache load failed for {self.uniprot_id}: {e}")
            return False


_CACHE: Dict[str, StructureFeatures] = {}


def get_structure_features(uniprot_id: str) -> StructureFeatures:
    """Process-level memoised loader."""
    sf = _CACHE.get(uniprot_id)
    if sf is None:
        sf = StructureFeatures(uniprot_id)
        _CACHE[uniprot_id] = sf
    return sf


if __name__ == "__main__":
    # Self-test on canonical GOF hotspots and a structural control.
    logging.basicConfig(level=logging.INFO)
    cases = [
        ("P42336", 1047, "PIK3CA H1047R — kinase activation hotspot"),
        ("P35222", 33,   "CTNNB1 S33 — phosphodegron site"),
        ("P05997", 1000, "COL5A2 — structural triple-helix control"),
        ("P04637", 175,  "TP53 R175 — DNA-binding (DN) hotspot"),
    ]
    for up, pos, label in cases:
        sf = get_structure_features(up)
        print(f"\n{label}  ({up}:{pos})  loaded={sf.ok}")
        if sf.has(pos):
            for k, v in sf.summary(pos).items():
                print(f"   {k}: {v}")

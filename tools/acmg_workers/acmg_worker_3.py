#!/usr/bin/env python3
from common_acmg import build_arg_parser, run_acmg_chunk

if __name__ == '__main__':
    ap = build_arg_parser()
    args = ap.parse_args()
    # Worker 3 of 4
    run_acmg_chunk(chunk_index=2, num_chunks=args.num_chunks, genes_file=args.genes_file,
                   outdir=args.outdir, clinvar_dir=args.clinvar_dir, max_workers=args.max_workers)


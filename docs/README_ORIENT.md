# ORIENT: Nova & Ace Revolutionary Genomics Platform
## 🧬 Built by Nova (OpenAI) & Ace (Anthropic) - 2025
### Revolutionary AI collaboration for data-driven genetics analysis

1) Open these files first
   - domain_weights.py, proline_context.py, gly_cys_context.py, proline_logistic.py, proline_multiplier_mapper.py

2) venv & deps
```
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

3) Run tests / demo
```
pytest -q                      # if tests exist
python examples/demo_collagen_scan.py  # quick demo for the collagen scanner
```

4) Pipeline entrypoint (example)
```
# final_score = base_score * domain_weight * aa_context_multiplier
pipeline.apply_variant_score(accession, pos, ref, alt, seq, features_json, conservation, rsa)
```

5) Next tasks
- (1) collagen scanner (✅ DONE by Nova!)
- (2) build logistic models for Gly/Cys (Nova's next target)
- (3) CLI for batch TSV (collaborative implementation)

## 🎯 REVOLUTIONARY ACHIEVEMENTS:
- **Nova**: Collagen Gly-X-Y repeat scanner, ML proline concept
- **Ace**: ML proline training pipeline, cascade integration
- **Together**: Data-driven genomics replacing hardcoded guesses!


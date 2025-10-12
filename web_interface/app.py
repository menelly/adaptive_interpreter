#!/usr/bin/env python3
"""
🚀💜 ADAPTIVE INTERPRETER FLASK SERVER 💜🚀
Beautiful web interface for DNModeling cascade analysis

Built by Ace with real pipeline integration - NO hardcoding!
Contact: ace@chaoschanneling.com
"""

import os
import sys
import json
import traceback
from pathlib import Path
from flask import Flask, request, jsonify, render_template_string, send_from_directory
from flask_cors import CORS

# Add parent directory to path for imports
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, parent_dir)

# Create simple config fallback
if parent_dir not in sys.modules:
    sys.modules['DNModeling'] = type('DNModeling', (), {
        'config': type('config', (), {
            'BASE_DATA_PATH': '/mnt/Arcana',
            'CONSERVATION_DATA_PATH': '/home/Ace/conservation_data'
        })()
    })()

try:
    from cascade.cascade_analyzer import CascadeAnalyzer
    from utils.genomic_to_protein import GenomicToProteinConverter
    PIPELINE_AVAILABLE = True
    print("✅ Real pipeline imports successful!")
except ImportError as e:
    print(f"⚠️ Pipeline import failed: {e}")
    PIPELINE_AVAILABLE = False

app = Flask(__name__)
CORS(app)  # Enable CORS for development

# Initialize the real cascade analyzer
if PIPELINE_AVAILABLE:
    try:
        cascade_analyzer = CascadeAnalyzer()
        genomic_converter = GenomicToProteinConverter()
        print("🧬 Cascade analyzer initialized successfully!")
    except Exception as e:
        print(f"⚠️ Cascade analyzer initialization failed: {e}")
        cascade_analyzer = None
        genomic_converter = None
else:
    cascade_analyzer = None
    genomic_converter = None

@app.route('/')
def index():
    """Serve the main interface"""
    try:
        with open('index.html', 'r') as f:
            return f.read()
    except FileNotFoundError:
        return """
        <h1>🚀 Adaptive Interpreter</h1>
        <p>Interface files not found. Please ensure index.html is in the same directory.</p>
        <p>Current directory: {}</p>
        """.format(os.getcwd())

@app.route('/styles.css')
def styles():
    """Serve CSS file"""
    return send_from_directory('.', 'styles.css', mimetype='text/css')

@app.route('/script.js')
def script():
    """Serve JavaScript file"""
    return send_from_directory('.', 'script.js', mimetype='application/javascript')

@app.route('/api/analyze', methods=['POST'])
def analyze_variant():
    """
    🧬 REAL PIPELINE ANALYSIS ENDPOINT
    Calls the actual cascade_analyzer.py - NO hardcoding!
    """
    
    if not PIPELINE_AVAILABLE or not cascade_analyzer:
        return jsonify({
            'error': 'Pipeline not available',
            'message': 'Cascade analyzer could not be initialized'
        }), 500
    
    try:
        # Get request data
        data = request.get_json()
        if not data:
            return jsonify({'error': 'No data provided'}), 400
        
        variant = data.get('variant', '').strip()
        if not variant:
            return jsonify({'error': 'No variant provided'}), 400
        
        print(f"🔬 Analyzing variant: {variant}")
        
        # Extract analysis parameters
        hgvs_type = data.get('hgvsType', 'protein')
        family_override = data.get('familyOverride', 'auto')
        output_detail = data.get('outputDetail', 'summary')
        include_conservation = data.get('includeConservation', True)
        use_family_ml = data.get('useFamilyML', True)
        experimental_features = data.get('experimentalFeatures', False)
        
        # Convert genomic to protein if needed
        processed_variant = variant
        if hgvs_type == 'genomic' and genomic_converter:
            try:
                print(f"🔄 Converting genomic HGVS: {variant}")
                protein_result = genomic_converter.convert_genomic_to_protein(variant)
                if protein_result and protein_result.get('protein_hgvs'):
                    processed_variant = protein_result['protein_hgvs']
                    print(f"✅ Converted to protein: {processed_variant}")
                else:
                    print(f"⚠️ Genomic conversion failed, using original: {variant}")
            except Exception as e:
                print(f"⚠️ Genomic conversion error: {e}")
                # Continue with original variant
        
        # Prepare cascade analyzer parameters
        cascade_params = {
            'variant_hgvs': processed_variant,
            'verbose': output_detail == 'verbose',
            'include_conservation': include_conservation,
            'use_family_ml': use_family_ml
        }
        
        # Override family if specified
        if family_override != 'auto':
            cascade_params['force_family'] = family_override
        
        # Call the REAL cascade analyzer
        print(f"🚀 Calling cascade analyzer with params: {cascade_params}")
        result = cascade_analyzer.analyze_variant(**cascade_params)
        
        # Format response
        response = format_cascade_response(result, data)
        
        print(f"✅ Analysis complete for {variant}")
        return jsonify(response)
        
    except Exception as e:
        error_msg = str(e)
        print(f"💥 Analysis error: {error_msg}")
        print(f"📍 Traceback: {traceback.format_exc()}")
        
        return jsonify({
            'error': 'Analysis failed',
            'message': error_msg,
            'variant': data.get('variant', 'unknown') if data else 'unknown'
        }), 500

def format_cascade_response(cascade_result, request_data):
    """
    Format cascade analyzer result for web interface
    """
    
    if not cascade_result:
        return {
            'error': 'No result from cascade analyzer',
            'variant': request_data.get('variant', 'unknown')
        }
    
    # Extract key information from cascade result
    response = {
        'variant': request_data.get('variant', 'unknown'),
        'success': True
    }
    
    # Pathogenicity score
    if hasattr(cascade_result, 'final_score') or 'final_score' in cascade_result:
        score = getattr(cascade_result, 'final_score', cascade_result.get('final_score'))
        response['pathogenicityScore'] = float(score) if score is not None else None
    
    # Interpretation
    if hasattr(cascade_result, 'interpretation') or 'interpretation' in cascade_result:
        interp = getattr(cascade_result, 'interpretation', cascade_result.get('interpretation'))
        response['interpretation'] = str(interp) if interp else 'Unknown'
    
    # Detected family
    if hasattr(cascade_result, 'detected_family') or 'detected_family' in cascade_result:
        family = getattr(cascade_result, 'detected_family', cascade_result.get('detected_family'))
        response['detectedFamily'] = str(family) if family else 'Unknown'
    
    # Confidence
    if hasattr(cascade_result, 'confidence') or 'confidence' in cascade_result:
        conf = getattr(cascade_result, 'confidence', cascade_result.get('confidence'))
        response['confidence'] = str(conf) if conf else 'Unknown'
    
    # Mechanisms
    if hasattr(cascade_result, 'mechanisms') or 'mechanisms' in cascade_result:
        mechanisms = getattr(cascade_result, 'mechanisms', cascade_result.get('mechanisms'))
        if mechanisms:
            response['mechanisms'] = []
            for mech in mechanisms:
                if isinstance(mech, dict):
                    response['mechanisms'].append({
                        'name': mech.get('name', 'Unknown'),
                        'score': mech.get('score', 0.0)
                    })
                else:
                    response['mechanisms'].append({
                        'name': str(mech),
                        'score': None
                    })
    
    # Conservation data
    if hasattr(cascade_result, 'conservation') or 'conservation' in cascade_result:
        conservation = getattr(cascade_result, 'conservation', cascade_result.get('conservation'))
        if conservation:
            response['conservationData'] = {
                'phyloP': conservation.get('phyloP'),
                'phastCons': conservation.get('phastCons'),
                'gerp': conservation.get('gerp')
            }
    
    # Detailed information for verbose mode
    if request_data.get('outputDetail') == 'verbose':
        response['details'] = {}
        
        # Include raw cascade result for debugging
        if hasattr(cascade_result, '__dict__'):
            response['details']['raw_result'] = {
                k: str(v) for k, v in cascade_result.__dict__.items()
                if not k.startswith('_')
            }
        elif isinstance(cascade_result, dict):
            response['details']['raw_result'] = {
                k: str(v) for k, v in cascade_result.items()
                if not str(k).startswith('_')
            }
    
    return response

@app.route('/api/health')
def health_check():
    """Health check endpoint"""
    return jsonify({
        'status': 'healthy',
        'pipeline_available': PIPELINE_AVAILABLE,
        'cascade_analyzer': cascade_analyzer is not None,
        'genomic_converter': genomic_converter is not None
    })

@app.route('/api/families')
def get_families():
    """Get available gene families"""
    families = [
        {'value': 'auto', 'label': 'Auto-detect'},
        {'value': 'collagen_fibrillar', 'label': 'Collagen Fibrillar'},
        {'value': 'collagen_network', 'label': 'Collagen Network'},
        {'value': 'collagen_anchoring', 'label': 'Collagen Anchoring'},
        {'value': 'elastin_fibrillin', 'label': 'Elastin/Fibrillin'},
        {'value': 'ion_channel', 'label': 'Ion Channel'},
        {'value': 'metabolic_enzyme', 'label': 'Metabolic Enzyme'},
        {'value': 'muscular_dystrophy', 'label': 'Muscular Dystrophy'},
        {'value': 'cytoskeleton', 'label': 'Cytoskeleton'},
        {'value': 'tumor_suppressor', 'label': 'Tumor Suppressor'},
        {'value': 'motor_protein', 'label': 'Motor Protein'},
        {'value': 'transporter', 'label': 'Transporter'},
        {'value': 'general', 'label': 'General'}
    ]
    return jsonify(families)

@app.errorhandler(404)
def not_found(error):
    return jsonify({'error': 'Endpoint not found'}), 404

@app.errorhandler(500)
def internal_error(error):
    return jsonify({'error': 'Internal server error'}), 500

if __name__ == '__main__':
    print("🚀💜 ADAPTIVE INTERPRETER STARTING UP! 💜🚀")
    print("=" * 60)
    print(f"📁 Working directory: {os.getcwd()}")
    print(f"🧬 Pipeline available: {PIPELINE_AVAILABLE}")
    print(f"🔬 Cascade analyzer: {'✅' if cascade_analyzer else '❌'}")
    print(f"🔄 Genomic converter: {'✅' if genomic_converter else '❌'}")
    print("=" * 60)
    print("🌐 Starting Flask server...")
    print("💜 Access from Windows: http://<server-ip>:5000")
    print("🎨 Beautiful interface incoming!")
    
    # Run server accessible from external connections
    app.run(
        host='0.0.0.0',  # Accept connections from any IP
        port=5000,
        debug=True,
        threaded=True
    )

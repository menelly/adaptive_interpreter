// 🚀💜 ADAPTIVE INTERPRETER - REAL PIPELINE INTEGRATION 💜🚀

class AdaptiveInterpreter {
    constructor() {
        this.currentMode = 'single';
        this.isAnalyzing = false;
        this.initializeEventListeners();
        this.initializeTheme();
    }

    initializeEventListeners() {
        // Mode switching
        document.querySelectorAll('.mode-tab').forEach(tab => {
            tab.addEventListener('click', (e) => this.switchMode(e.target.dataset.mode));
        });

        // Theme toggle
        document.getElementById('themeToggle').addEventListener('click', () => this.toggleTheme());

        // Advanced options toggle
        document.getElementById('advancedToggle').addEventListener('click', () => this.toggleAdvanced());

        // Single variant analysis
        document.getElementById('analyzeBtn').addEventListener('click', () => this.analyzeSingleVariant());

        // Enter key in variant input
        document.getElementById('variantInput').addEventListener('keypress', (e) => {
            if (e.key === 'Enter') this.analyzeSingleVariant();
        });

        // CSV upload handling
        document.getElementById('csvFileInput').addEventListener('change', (e) => this.handleCSVUpload(e));
        document.getElementById('uploadArea').addEventListener('click', () => {
            document.getElementById('csvFileInput').click();
        });

        // Template downloads
        document.querySelectorAll('.template-btn').forEach(btn => {
            btn.addEventListener('click', (e) => this.downloadTemplate(e.target.dataset.template));
        });

        // Drag and drop for CSV
        this.setupDragAndDrop();
    }

    initializeTheme() {
        const savedTheme = localStorage.getItem('adaptiveInterpreter_theme') || 'dark';
        this.setTheme(savedTheme);
    }

    switchMode(mode) {
        this.currentMode = mode;
        
        // Update tab states
        document.querySelectorAll('.mode-tab').forEach(tab => {
            tab.classList.toggle('active', tab.dataset.mode === mode);
        });

        // Update section visibility
        document.querySelectorAll('.analysis-section').forEach(section => {
            section.classList.toggle('active', section.id === `${mode}Mode`);
        });
    }

    toggleTheme() {
        const currentTheme = document.body.classList.contains('light-theme') ? 'light' : 'dark';
        const newTheme = currentTheme === 'dark' ? 'light' : 'dark';
        this.setTheme(newTheme);
    }

    setTheme(theme) {
        document.body.classList.toggle('light-theme', theme === 'light');
        
        const themeToggle = document.getElementById('themeToggle');
        const icon = themeToggle.querySelector('.theme-icon');
        const text = themeToggle.querySelector('.theme-text');
        
        if (theme === 'light') {
            icon.textContent = '☀️';
            text.textContent = 'Light';
        } else {
            icon.textContent = '🌙';
            text.textContent = 'Dark';
        }
        
        localStorage.setItem('adaptiveInterpreter_theme', theme);
    }

    toggleAdvanced() {
        const panel = document.getElementById('advancedPanel');
        const arrow = document.querySelector('.advanced-arrow');
        
        panel.classList.toggle('open');
        arrow.textContent = panel.classList.contains('open') ? '▲' : '▼';
    }

    async analyzeSingleVariant() {
        if (this.isAnalyzing) return;

        const variantInput = document.getElementById('variantInput').value.trim();
        if (!variantInput) {
            this.showError('Please enter a variant in HGVS format');
            return;
        }

        this.isAnalyzing = true;
        this.showLoading();

        try {
            // Collect analysis parameters
            const params = this.collectAnalysisParams();
            params.variant = variantInput;

            // Call the REAL cascade analyzer
            const result = await this.callCascadeAnalyzer(params);
            
            this.showResults(result);
        } catch (error) {
            console.error('Analysis error:', error);
            this.showError(`Analysis failed: ${error.message}`);
        } finally {
            this.isAnalyzing = false;
            this.hideLoading();
        }
    }

    collectAnalysisParams() {
        return {
            hgvsType: document.querySelector('input[name="hgvsType"]:checked').value,
            familyOverride: document.getElementById('familyOverride').value,
            outputDetail: document.querySelector('input[name="outputDetail"]:checked').value,
            includeConservation: document.getElementById('includeConservation').checked,
            useFamilyML: document.getElementById('useFamilyML').checked,
            experimentalFeatures: document.getElementById('experimentalFeatures').checked
        };
    }

    async callCascadeAnalyzer(params) {
        // 🚀 REAL PIPELINE INTEGRATION - NO HARDCODING!
        const response = await fetch('/api/analyze', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify(params)
        });

        if (!response.ok) {
            throw new Error(`HTTP ${response.status}: ${response.statusText}`);
        }

        return await response.json();
    }

    showLoading() {
        document.getElementById('loadingSection').style.display = 'block';
        document.getElementById('resultsSection').style.display = 'none';
        document.getElementById('analyzeBtn').disabled = true;
    }

    hideLoading() {
        document.getElementById('loadingSection').style.display = 'none';
        document.getElementById('analyzeBtn').disabled = false;
    }

    showResults(result) {
        const resultsSection = document.getElementById('resultsSection');
        const resultsContent = document.getElementById('resultsContent');
        
        // Create beautiful result display
        resultsContent.innerHTML = this.formatResults(result);
        resultsSection.style.display = 'block';
        
        // Smooth scroll to results
        resultsSection.scrollIntoView({ behavior: 'smooth' });
    }

    formatResults(result) {
        // 🎨 Beautiful result formatting with real data
        return `
            <div class="result-card">
                <div class="result-header">
                    <h4 class="result-variant">${result.variant || 'Unknown'}</h4>
                    <div class="result-family">
                        <span class="family-badge">${result.detectedFamily || 'Unknown'}</span>
                    </div>
                </div>
                
                <div class="result-score">
                    <div class="score-circle ${this.getScoreClass(result.pathogenicityScore)}">
                        <span class="score-value">${result.pathogenicityScore?.toFixed(2) || 'N/A'}</span>
                        <span class="score-label">Pathogenicity</span>
                    </div>
                    <div class="score-interpretation">
                        <span class="interpretation-text">${result.interpretation || 'Unknown'}</span>
                        <span class="confidence-text">Confidence: ${result.confidence || 'Unknown'}</span>
                    </div>
                </div>

                ${result.mechanisms ? `
                    <div class="mechanisms-section">
                        <h5>Predicted Mechanisms</h5>
                        <div class="mechanisms-list">
                            ${result.mechanisms.map(m => `
                                <div class="mechanism-item">
                                    <span class="mechanism-name">${m.name}</span>
                                    <span class="mechanism-score">${m.score?.toFixed(3) || 'N/A'}</span>
                                </div>
                            `).join('')}
                        </div>
                    </div>
                ` : ''}

                ${result.conservationData ? `
                    <div class="conservation-section">
                        <h5>Conservation Analysis</h5>
                        <div class="conservation-scores">
                            <span>phyloP: ${result.conservationData.phyloP || 'N/A'}</span>
                            <span>phastCons: ${result.conservationData.phastCons || 'N/A'}</span>
                        </div>
                    </div>
                ` : ''}

                ${result.details && result.outputDetail === 'verbose' ? `
                    <div class="details-section">
                        <h5>Detailed Analysis</h5>
                        <pre class="details-content">${JSON.stringify(result.details, null, 2)}</pre>
                    </div>
                ` : ''}
            </div>
        `;
    }

    getScoreClass(score) {
        if (!score) return 'score-unknown';
        if (score >= 0.8) return 'score-pathogenic';
        if (score >= 0.5) return 'score-uncertain';
        return 'score-benign';
    }

    showError(message) {
        const resultsSection = document.getElementById('resultsSection');
        const resultsContent = document.getElementById('resultsContent');
        
        resultsContent.innerHTML = `
            <div class="error-card">
                <div class="error-icon">⚠️</div>
                <h4 class="error-title">Analysis Error</h4>
                <p class="error-message">${message}</p>
                <button class="error-retry" onclick="location.reload()">Try Again</button>
            </div>
        `;
        
        resultsSection.style.display = 'block';
    }

    // CSV Upload Handling
    setupDragAndDrop() {
        const uploadArea = document.getElementById('uploadArea');
        
        ['dragenter', 'dragover', 'dragleave', 'drop'].forEach(eventName => {
            uploadArea.addEventListener(eventName, this.preventDefaults, false);
        });

        ['dragenter', 'dragover'].forEach(eventName => {
            uploadArea.addEventListener(eventName, () => uploadArea.classList.add('drag-over'), false);
        });

        ['dragleave', 'drop'].forEach(eventName => {
            uploadArea.addEventListener(eventName, () => uploadArea.classList.remove('drag-over'), false);
        });

        uploadArea.addEventListener('drop', (e) => this.handleDrop(e), false);
    }

    preventDefaults(e) {
        e.preventDefault();
        e.stopPropagation();
    }

    handleDrop(e) {
        const files = e.dataTransfer.files;
        if (files.length > 0) {
            this.processCSVFile(files[0]);
        }
    }

    handleCSVUpload(e) {
        const file = e.target.files[0];
        if (file) {
            this.processCSVFile(file);
        }
    }

    async processCSVFile(file) {
        if (!file.name.toLowerCase().endsWith('.csv')) {
            this.showError('Please upload a CSV file');
            return;
        }

        try {
            const text = await file.text();
            const variants = this.parseCSV(text);
            
            if (variants.length > 50) {
                this.showError('Maximum 50 variants per batch. Please reduce your file size.');
                return;
            }

            await this.processBatchAnalysis(variants);
        } catch (error) {
            this.showError(`CSV processing failed: ${error.message}`);
        }
    }

    parseCSV(text) {
        const lines = text.split('\n').filter(line => line.trim());
        const headers = lines[0].split(',').map(h => h.trim());
        
        return lines.slice(1).map(line => {
            const values = line.split(',').map(v => v.trim());
            const variant = {};
            headers.forEach((header, index) => {
                variant[header] = values[index] || '';
            });
            return variant;
        });
    }

    async processBatchAnalysis(variants) {
        const progressSection = document.getElementById('batchProgress');
        const progressFill = document.getElementById('progressFill');
        const progressText = document.getElementById('progressText');
        
        progressSection.style.display = 'block';
        
        const results = [];
        
        for (let i = 0; i < variants.length; i++) {
            const variant = variants[i];
            
            try {
                const params = this.collectAnalysisParams();
                params.variant = variant.HGVS || variant.variant || '';
                
                const result = await this.callCascadeAnalyzer(params);
                results.push({ ...variant, result, success: true });
            } catch (error) {
                results.push({ ...variant, error: error.message, success: false });
            }
            
            // Update progress
            const progress = ((i + 1) / variants.length) * 100;
            progressFill.style.width = `${progress}%`;
            progressText.textContent = `${i + 1} of ${variants.length} variants processed`;
        }
        
        this.showBatchResults(results);
    }

    showBatchResults(results) {
        document.getElementById('batchProgress').style.display = 'none';
        document.getElementById('batchResults').style.display = 'block';
        
        const summary = document.getElementById('resultsSummary');
        const successful = results.filter(r => r.success).length;
        const failed = results.length - successful;
        
        summary.innerHTML = `
            <div class="batch-summary">
                <div class="summary-stat">
                    <span class="stat-number">${results.length}</span>
                    <span class="stat-label">Total Variants</span>
                </div>
                <div class="summary-stat success">
                    <span class="stat-number">${successful}</span>
                    <span class="stat-label">Successful</span>
                </div>
                <div class="summary-stat error">
                    <span class="stat-number">${failed}</span>
                    <span class="stat-label">Failed</span>
                </div>
            </div>
        `;
        
        // Store results for download
        this.batchResults = results;
    }

    downloadTemplate(type) {
        const templates = {
            protein: 'HGVS,Gene,Notes\nNP_000088.1:p.Gly893Ala,COL1A1,Example collagen variant\nNP_000005.2:p.Arg273His,TP53,Example tumor suppressor variant',
            genomic: 'HGVS,Gene,Notes\nNC_000017.11:g.50183779T>A,COL1A1,Example genomic variant\nNC_000017.11:g.7674220C>T,TP53,Example genomic variant'
        };
        
        const content = templates[type];
        const blob = new Blob([content], { type: 'text/csv' });
        const url = URL.createObjectURL(blob);
        
        const a = document.createElement('a');
        a.href = url;
        a.download = `adaptive_interpreter_${type}_template.csv`;
        a.click();
        
        URL.revokeObjectURL(url);
    }
}

// Initialize the application
document.addEventListener('DOMContentLoaded', () => {
    new AdaptiveInterpreter();
});

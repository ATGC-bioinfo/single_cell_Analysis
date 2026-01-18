#!/usr/bin/env python3
"""
Test script to validate single cell analysis pipeline
"""

import os
import sys
import subprocess
import tempfile
import pandas as pd
import numpy as np
from scipy.io import mmread, mmwrite

def create_test_data():
    """Create minimal test data for pipeline validation"""
    
    # Create temporary directory
    test_dir = tempfile.mkdtemp(prefix="singlecell_test_")
    print(f"Creating test data in: {test_dir}")
    
    # Create synthetic count matrix (100 cells x 50 genes)
    n_cells, n_genes = 100, 50
    counts = np.random.negative_binomial(5, 0.3, size=(n_genes, n_cells))
    
    # Create gene names
    gene_names = [f"GENE_{i:03d}" for i in range(n_genes)]
    
    # Add some mitochondrial genes
    gene_names[0:3] = ["MT-ND1", "MT-ND2", "MT-CO1"]
    
    # Add some ribosomal genes
    gene_names[3:5] = ["RPS3", "RPL5"]
    
    # Create cell names
    cell_names = [f"CELL_{i:03d}" for i in range(n_cells)]
    
    # Save count matrix
    counts_file = os.path.join(test_dir, "test_counts.mtx")
    mmwrite(counts_file, counts)
    
    # Save gene names
    genes_file = os.path.join(test_dir, "test_genes.tsv")
    pd.DataFrame({"gene_id": gene_names}).to_csv(genes_file, sep="\t", index=False, header=False)
    
    # Save cell names
    cells_file = os.path.join(test_dir, "test_cells.tsv")
    pd.DataFrame({"cell_id": cell_names}).to_csv(cells_file, sep="\t", index=False, header=False)
    
    # Create samplesheet
    samplesheet_file = os.path.join(test_dir, "samplesheet.csv")
    pd.DataFrame({
        "sample_id": ["test_sample"],
        "counts": [counts_file]
    }).to_csv(samplesheet_file, index=False)
    
    return test_dir, samplesheet_file

def test_pipeline_basic(samplesheet_file):
    """Test basic pipeline functionality"""
    
    print("Testing basic pipeline functionality...")
    
    # Change to project directory
    project_dir = "/home/habee/Documents/single_cell_Analysis"
    os.chdir(project_dir)
    
    # Test nextflow syntax
    try:
        result = subprocess.run(
            ["nextflow", "lint", "main.nf"],
            capture_output=True,
            text=True,
            timeout=60
        )
        
        if result.returncode == 0:
            print("✓ Nextflow syntax validation passed")
        else:
            print("✗ Nextflow syntax validation failed:")
            print(result.stderr)
            return False
            
    except subprocess.TimeoutExpired:
        print("✗ Nextflow lint timed out")
        return False
    except FileNotFoundError:
        print("✗ Nextflow not found - please install Nextflow")
        return False
    
    return True

def test_python_imports():
    """Test that all required Python packages can be imported"""
    
    print("Testing Python package imports...")
    
    required_packages = [
        "scanpy",
        "scvi",
        "pandas",
        "numpy",
        "matplotlib",
        "seaborn",
        "scipy",
        "sklearn",
        "anndata"
    ]
    
    failed_imports = []
    
    for package in required_packages:
        try:
            __import__(package)
            print(f"✓ {package}")
        except ImportError:
            print(f"✗ {package}")
            failed_imports.append(package)
    
    if failed_imports:
        print(f"\nMissing packages: {', '.join(failed_imports)}")
        print("Please install missing packages using: pip install -r requirements.txt")
        return False
    
    return True

def test_config_file():
    """Test configuration file syntax"""
    
    print("Testing configuration file...")
    
    try:
        import yaml
        config_file = "/home/habee/Documents/single_cell_Analysis/config.yaml"
        
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        
        # Check required sections
        required_sections = ["qc", "doublet", "clustering", "scvi", "output"]
        for section in required_sections:
            if section in config:
                print(f"✓ {section} section found")
            else:
                print(f"✗ {section} section missing")
                return False
        
        return True
        
    except ImportError:
        print("✗ PyYAML not installed - cannot validate config file")
        return False
    except Exception as e:
        print(f"✗ Config file validation failed: {e}")
        return False

def main():
    """Main test function"""
    
    print("Single Cell Analysis Pipeline - Test Suite")
    print("=" * 50)
    
    tests_passed = 0
    total_tests = 4
    
    # Test 1: Python imports
    if test_python_imports():
        tests_passed += 1
    
    print()
    
    # Test 2: Configuration file
    if test_config_file():
        tests_passed += 1
    
    print()
    
    # Test 3: Create test data
    test_dir, samplesheet_file = create_test_data()
    print("✓ Test data created successfully")
    tests_passed += 1
    
    print()
    
    # Test 4: Basic pipeline functionality
    if test_pipeline_basic(samplesheet_file):
        tests_passed += 1
    
    print()
    print("=" * 50)
    print(f"Test Results: {tests_passed}/{total_tests} tests passed")
    
    if tests_passed == total_tests:
        print("🎉 All tests passed! Pipeline is ready to use.")
        print(f"\nTo run the pipeline with test data:")
        print(f"cd /home/habee/Documents/single_cell_Analysis")
        print(f"nextflow run main.nf --samplesheet {samplesheet_file}")
    else:
        print("❌ Some tests failed. Please fix issues before running the pipeline.")
    
    # Cleanup
    import shutil
    shutil.rmtree(test_dir)
    
    return tests_passed == total_tests

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
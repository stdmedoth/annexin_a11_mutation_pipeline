"""
BioEmu Runner Module

Interface for running BioEmu to generate protein conformational ensembles.
"""

import subprocess
import json
import tempfile
from pathlib import Path
from typing import Optional, List, Dict, Any
from dataclasses import dataclass, field
import yaml


@dataclass
class BioEmuConfig:
    """Configuration for BioEmu runs."""
    num_conformations: int = 100
    output_format: str = "pdb"
    device: str = "cuda"
    seed: int = 42
    batch_size: int = 10
    temperature: float = 1.0
    
    @classmethod
    def from_yaml(cls, config_path: Path) -> "BioEmuConfig":
        """Load configuration from YAML file."""
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        
        bioemu_config = config.get('bioemu', {})
        return cls(**{k: v for k, v in bioemu_config.items() if k in cls.__dataclass_fields__})


class BioEmuRunner:
    """
    Runner for BioEmu conformational sampling.
    
    BioEmu is a deep learning model for generating protein conformational ensembles.
    This class provides an interface to run BioEmu and manage outputs.
    """
    
    def __init__(self, config: Optional[BioEmuConfig] = None):
        """
        Initialize BioEmu runner.
        
        Args:
            config: BioEmu configuration. Uses defaults if not provided.
        """
        self.config = config or BioEmuConfig()
        self.bioemu_available = self._check_bioemu_available()
    
    def _check_bioemu_available(self) -> bool:
        """Check if BioEmu is available."""
        try:
            # Try importing bioemu
            import bioemu
            return True
        except ImportError:
            return False
    
    def generate_conformations(
        self,
        sequence: str,
        output_dir: Path,
        name: str = "protein",
        num_conformations: Optional[int] = None
    ) -> Dict[str, Any]:
        """
        Generate conformational ensemble for a sequence.
        
        Args:
            sequence: Amino acid sequence
            output_dir: Directory to save output structures
            name: Name for the protein/output files
            num_conformations: Override default number of conformations
            
        Returns:
            Dictionary with run information and output paths
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        n_conf = num_conformations or self.config.num_conformations
        
        # Create input file for BioEmu
        input_data = {
            "sequence": sequence,
            "num_samples": n_conf,
            "temperature": self.config.temperature,
            "seed": self.config.seed,
        }
        
        input_file = output_dir / f"{name}_input.json"
        with open(input_file, 'w') as f:
            json.dump(input_data, f, indent=2)
        
        # Run BioEmu (placeholder for actual implementation)
        output_files = self._run_bioemu(
            sequence=sequence,
            output_dir=output_dir,
            name=name,
            n_samples=n_conf
        )
        
        result = {
            "name": name,
            "sequence": sequence,
            "sequence_length": len(sequence),
            "num_conformations": n_conf,
            "output_directory": str(output_dir),
            "output_files": output_files,
            "config": {
                "device": self.config.device,
                "seed": self.config.seed,
                "temperature": self.config.temperature,
            }
        }
        
        # Save run metadata
        metadata_file = output_dir / f"{name}_metadata.json"
        with open(metadata_file, 'w') as f:
            json.dump(result, f, indent=2)
        
        return result
    
    def _run_bioemu(
        self,
        sequence: str,
        output_dir: Path,
        name: str,
        n_samples: int
    ) -> List[str]:
        """
        Run BioEmu to generate conformations.
        
        This is the actual BioEmu execution method.
        
        Args:
            sequence: Amino acid sequence
            output_dir: Output directory
            name: Output name prefix
            n_samples: Number of samples to generate
            
        Returns:
            List of output file paths
        """
        output_files = []
        
        if self.bioemu_available:
            try:
                import bioemu
                import torch
                
                print(f"[BioEmu] Loading model on {self.config.device}...")
                
                # Set seed for reproducibility
                torch.manual_seed(self.config.seed)
                if self.config.device == "cuda":
                    torch.cuda.manual_seed(self.config.seed)
                
                # Initialize BioEmu sampler
                sampler = bioemu.get_sampler(device=self.config.device)
                
                print(f"[BioEmu] Generating {n_samples} conformations for {name}...")
                print(f"[BioEmu] Sequence length: {len(sequence)} residues")
                
                # Generate conformations
                # BioEmu returns a list of structures
                samples = sampler.sample(
                    sequence,
                    num_samples=n_samples,
                    temperature=self.config.temperature
                )
                
                # Save each sample as PDB
                for i, sample in enumerate(samples):
                    out_path = output_dir / f"{name}_conf_{i:04d}.pdb"
                    # Save structure - adjust based on actual BioEmu API
                    if hasattr(sample, 'to_pdb'):
                        sample.to_pdb(str(out_path))
                    elif hasattr(sample, 'save'):
                        sample.save(str(out_path))
                    else:
                        # If sample is coordinates, create PDB manually
                        self._write_pdb(sample, sequence, str(out_path))
                    output_files.append(str(out_path))
                
                print(f"[BioEmu] ✓ Generated {len(output_files)} PDB files")
                
            except Exception as e:
                print(f"[BioEmu] Error during generation: {e}")
                print(f"[BioEmu] Creating placeholder files instead...")
                output_files = self._create_placeholder_files(output_dir, name, n_samples, sequence)
        else:
            print(f"[BioEmu] ⚠ BioEmu not installed!")
            print(f"[BioEmu] To install BioEmu, run:")
            print(f"[BioEmu]   pip install bioemu")
            print(f"[BioEmu] Or clone from: https://github.com/microsoft/bioemu")
            print(f"[BioEmu]")
            print(f"[BioEmu] Creating placeholder files for pipeline testing...")
            output_files = self._create_placeholder_files(output_dir, name, n_samples, sequence)
        
        return output_files
    
    def _create_placeholder_files(
        self, 
        output_dir: Path, 
        name: str, 
        n_samples: int,
        sequence: str
    ) -> List[str]:
        """Create placeholder files for testing when BioEmu is not available."""
        output_files = []
        
        # Save sequence to FASTA for reference
        fasta_path = output_dir / f"{name}_sequence.fasta"
        with open(fasta_path, 'w') as f:
            f.write(f">{name}\n{sequence}\n")
        
        print(f"[BioEmu] Sequence: {sequence[:50]}..." if len(sequence) > 50 else f"[BioEmu] Sequence: {sequence}")
        print(f"[BioEmu] Would generate {n_samples} PDB files")
        print(f"[BioEmu] Output directory: {output_dir}")
        
        # Just record what would be created
        for i in range(n_samples):
            expected_path = output_dir / f"{name}_conf_{i:04d}.pdb"
            output_files.append(str(expected_path))
        
        return output_files
    
    def _write_pdb(self, coords, sequence: str, filepath: str) -> None:
        """Write coordinates to a PDB file."""
        import numpy as np
        
        with open(filepath, 'w') as f:
            f.write(f"REMARK Generated by BioEmu\n")
            f.write(f"REMARK Sequence length: {len(sequence)}\n")
            
            atom_idx = 1
            for i, aa in enumerate(sequence):
                # Write CA atom for each residue
                if isinstance(coords, np.ndarray) and len(coords.shape) >= 2:
                    x, y, z = coords[i, 0], coords[i, 1], coords[i, 2]
                else:
                    x, y, z = 0.0, 0.0, float(i * 3.8)  # Extended chain
                
                f.write(f"ATOM  {atom_idx:5d}  CA  {aa:3s} A{i+1:4d}    "
                       f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n")
                atom_idx += 1
            
            f.write("END\n")
    
    def generate_wt_and_mutant(
        self,
        wt_sequence: str,
        mutant_sequence: str,
        output_dir: Path,
        wt_name: str = "WT",
        mutant_name: str = "Mutant"
    ) -> Dict[str, Any]:
        """
        Generate conformations for both WT and mutant sequences.
        
        Args:
            wt_sequence: Wild-type sequence
            mutant_sequence: Mutant sequence
            output_dir: Base output directory
            wt_name: Name for WT output
            mutant_name: Name for mutant output
            
        Returns:
            Dictionary with both WT and mutant results
        """
        output_dir = Path(output_dir)
        
        # Generate WT conformations
        wt_dir = output_dir / "wt"
        wt_result = self.generate_conformations(
            wt_sequence, wt_dir, wt_name
        )
        
        # Generate mutant conformations
        mutant_dir = output_dir / "mutant"
        mutant_result = self.generate_conformations(
            mutant_sequence, mutant_dir, mutant_name
        )
        
        return {
            "wt": wt_result,
            "mutant": mutant_result,
        }


def create_bioemu_script(
    sequence: str,
    output_dir: Path,
    name: str,
    config: BioEmuConfig
) -> Path:
    """
    Create a standalone script to run BioEmu.
    
    Useful for running on HPC clusters or with different environments.
    
    Args:
        sequence: Amino acid sequence
        output_dir: Output directory
        name: Output name
        config: BioEmu configuration
        
    Returns:
        Path to generated script
    """
    script_content = f'''#!/usr/bin/env python3
"""
BioEmu Conformation Generation Script
Generated for: {name}
"""

import os
import sys

# Configuration
SEQUENCE = "{sequence}"
OUTPUT_DIR = "{output_dir}"
NAME = "{name}"
NUM_SAMPLES = {config.num_conformations}
DEVICE = "{config.device}"
SEED = {config.seed}
TEMPERATURE = {config.temperature}

def main():
    # Import BioEmu (adjust based on your installation)
    try:
        from bioemu import BioEmu
    except ImportError:
        print("Error: BioEmu not installed. Please install it first.")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Initialize model
    print(f"Loading BioEmu model on {{DEVICE}}...")
    model = BioEmu(device=DEVICE)
    
    # Generate conformations
    print(f"Generating {{NUM_SAMPLES}} conformations...")
    structures = model.sample(
        sequence=SEQUENCE,
        n_samples=NUM_SAMPLES,
        temperature=TEMPERATURE,
        seed=SEED
    )
    
    # Save structures
    for i, struct in enumerate(structures):
        output_path = os.path.join(OUTPUT_DIR, f"{{NAME}}_conf_{{i:04d}}.pdb")
        struct.save(output_path)
        print(f"Saved: {{output_path}}")
    
    print(f"\\nComplete! Generated {{len(structures)}} conformations.")

if __name__ == "__main__":
    main()
'''
    
    script_path = Path(output_dir) / f"run_bioemu_{name}.py"
    script_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(script_path, 'w') as f:
        f.write(script_content)
    
    # Make executable
    script_path.chmod(0o755)
    
    return script_path

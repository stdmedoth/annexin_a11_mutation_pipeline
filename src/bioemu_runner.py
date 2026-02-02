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
        self._check_bioemu_available()
    
    def _check_bioemu_available(self) -> bool:
        """Check if BioEmu is available."""
        # This is a placeholder - actual implementation depends on BioEmu installation
        # BioEmu might be a Python package, Docker container, or API
        self.bioemu_available = True
        return self.bioemu_available
    
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
        Modify based on your BioEmu installation method.
        
        Args:
            sequence: Amino acid sequence
            output_dir: Output directory
            name: Output name prefix
            n_samples: Number of samples to generate
            
        Returns:
            List of output file paths
        """
        output_files = []
        
        try:
            # Option 1: If BioEmu is installed as a Python package
            # from bioemu import BioEmu
            # model = BioEmu(device=self.config.device)
            # structures = model.sample(sequence, n_samples=n_samples)
            # for i, struct in enumerate(structures):
            #     out_path = output_dir / f"{name}_conf_{i:04d}.pdb"
            #     struct.save(out_path)
            #     output_files.append(str(out_path))
            
            # Option 2: If BioEmu is a command-line tool
            # cmd = [
            #     "bioemu",
            #     "--sequence", sequence,
            #     "--n_samples", str(n_samples),
            #     "--output_dir", str(output_dir),
            #     "--output_prefix", name,
            #     "--device", self.config.device,
            #     "--seed", str(self.config.seed),
            # ]
            # subprocess.run(cmd, check=True)
            
            # Option 3: If BioEmu runs via Docker
            # cmd = [
            #     "docker", "run", "--gpus", "all",
            #     "-v", f"{output_dir}:/output",
            #     "bioemu:latest",
            #     "--sequence", sequence,
            #     "--n_samples", str(n_samples),
            # ]
            # subprocess.run(cmd, check=True)
            
            # Placeholder: Create a template showing expected output structure
            print(f"[BioEmu] Generating {n_samples} conformations for {name}...")
            print(f"[BioEmu] Sequence length: {len(sequence)} residues")
            print(f"[BioEmu] Output directory: {output_dir}")
            
            # In actual implementation, BioEmu generates PDB files
            # For now, log what would happen
            for i in range(n_samples):
                expected_path = output_dir / f"{name}_conf_{i:04d}.pdb"
                output_files.append(str(expected_path))
            
            print(f"[BioEmu] Would generate {len(output_files)} PDB files")
            
            # Save sequence to FASTA for reference
            fasta_path = output_dir / f"{name}_sequence.fasta"
            with open(fasta_path, 'w') as f:
                f.write(f">{name}\n{sequence}\n")
            
        except Exception as e:
            print(f"[BioEmu] Error: {e}")
            raise RuntimeError(f"BioEmu execution failed: {e}")
        
        return output_files
    
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

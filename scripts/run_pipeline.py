#!/usr/bin/env python3
"""
Mutation Conformation Pipeline - Main Script

Interactive pipeline for generating and analyzing protein conformations
for mutations using BioEmu.

Usage:
    python scripts/run_pipeline.py
    python scripts/run_pipeline.py --config config/settings.yaml
"""

import sys
import os
from pathlib import Path
from datetime import datetime
from typing import Optional, List
import json
import yaml

# Add src to path
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.sequence_handler import SequenceHandler, validate_sequence_input, get_sequence_stats
from src.mutation_handler import MutationHandler, Mutation, create_mutation_directory
from src.bioemu_runner import BioEmuRunner, BioEmuConfig, create_bioemu_script
from src.alignment import StructureAligner, ConformationManager, AlignmentConfig


# Terminal colors
class Colors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    END = '\033[0m'
    BOLD = '\033[1m'


def print_header(text: str) -> None:
    """Print a formatted header."""
    print(f"\n{Colors.HEADER}{Colors.BOLD}{'=' * 60}{Colors.END}")
    print(f"{Colors.HEADER}{Colors.BOLD}{text:^60}{Colors.END}")
    print(f"{Colors.HEADER}{Colors.BOLD}{'=' * 60}{Colors.END}\n")


def print_success(text: str) -> None:
    """Print success message."""
    print(f"{Colors.GREEN}✓ {text}{Colors.END}")


def print_error(text: str) -> None:
    """Print error message."""
    print(f"{Colors.FAIL}✗ {text}{Colors.END}")


def print_info(text: str) -> None:
    """Print info message."""
    print(f"{Colors.CYAN}ℹ {text}{Colors.END}")


def print_warning(text: str) -> None:
    """Print warning message."""
    print(f"{Colors.WARNING}⚠ {text}{Colors.END}")


def load_config(config_path: Optional[Path] = None) -> dict:
    """Load configuration from YAML file."""
    if config_path is None:
        config_path = PROJECT_ROOT / "config" / "settings.yaml"
    
    if config_path.exists():
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    return {}


def get_sequence_input() -> SequenceHandler:
    """
    Interactive input for wild-type sequence.
    
    Returns:
        SequenceHandler with validated sequence
    """
    print_header("WILD-TYPE SEQUENCE INPUT")
    
    print("Please provide the wild-type protein sequence.")
    print("You can:")
    print("  1. Paste the sequence directly (one-letter amino acid code)")
    print("  2. Enter a file path to a FASTA file")
    print()
    
    while True:
        user_input = input(f"{Colors.BOLD}Enter sequence or FASTA path: {Colors.END}").strip()
        
        if not user_input:
            print_error("Input cannot be empty. Please try again.")
            continue
        
        # Check if it's a file path
        input_path = Path(user_input)
        if input_path.exists() and input_path.suffix.lower() in ['.fasta', '.fa', '.faa']:
            try:
                handler = SequenceHandler.from_fasta(input_path)
                print_success(f"Loaded sequence '{handler.name}' from file")
                print_info(f"Sequence length: {len(handler)} amino acids")
                return handler
            except Exception as e:
                print_error(f"Failed to load FASTA file: {e}")
                continue
        
        # Try to parse as direct sequence
        is_valid, error = validate_sequence_input(user_input)
        
        if is_valid:
            # Ask for a name
            name = input(f"{Colors.BOLD}Enter a name for this protein [protein]: {Colors.END}").strip()
            if not name:
                name = "protein"
            
            handler = SequenceHandler(user_input, name)
            print_success(f"Sequence validated successfully!")
            print_info(f"Name: {handler.name}")
            print_info(f"Length: {len(handler)} amino acids")
            
            # Show statistics
            stats = get_sequence_stats(str(handler))
            print_info(f"Charged residues: {stats['charged_residues']}")
            print_info(f"Hydrophobic residues: {stats['hydrophobic_residues']}")
            
            return handler
        else:
            print_error(error)
            print("Please try again.")


def get_mutation_input(wt_handler: SequenceHandler) -> Mutation:
    """
    Interactive input for mutation.
    
    Args:
        wt_handler: Wild-type sequence handler
        
    Returns:
        Validated Mutation object
    """
    print_header("MUTATION INPUT")
    
    print("Please specify the mutation to study.")
    print("Format: [OriginalAA][Position][MutantAA]")
    print("Examples:")
    print("  A123G  - Alanine at position 123 to Glycine")
    print("  V456L  - Valine at position 456 to Leucine")
    print()
    print(f"Sequence length: {len(wt_handler)} residues")
    print()
    
    mutation_handler = MutationHandler(wt_handler)
    
    while True:
        mutation_str = input(f"{Colors.BOLD}Enter mutation: {Colors.END}").strip().upper()
        
        if not mutation_str:
            print_error("Mutation cannot be empty. Please try again.")
            continue
        
        try:
            mutation = mutation_handler.add_mutation_from_string(mutation_str)
            
            print_success(f"Mutation validated: {mutation}")
            print_info(f"Original: {mutation.original} at position {mutation.position}")
            print_info(f"Mutant: {mutation.mutant}")
            
            # Show local sequence context
            start = max(0, mutation.position - 6)
            end = min(len(wt_handler), mutation.position + 5)
            context = str(wt_handler)[start:end]
            relative_pos = mutation.position - start - 1
            
            context_display = context[:relative_pos] + f"[{context[relative_pos]}]" + context[relative_pos + 1:]
            print_info(f"Local context: ...{context_display}...")
            
            return mutation
            
        except ValueError as e:
            print_error(str(e))
            print("Please try again.")


def get_additional_mutations(
    wt_handler: SequenceHandler,
    existing_mutations: List[Mutation]
) -> List[Mutation]:
    """
    Ask user if they want to add more mutations.
    
    Args:
        wt_handler: Wild-type sequence handler
        existing_mutations: List of already defined mutations
        
    Returns:
        Updated list of mutations
    """
    mutations = existing_mutations.copy()
    
    while True:
        add_more = input(f"\n{Colors.BOLD}Add another mutation? (y/N): {Colors.END}").strip().lower()
        
        if add_more not in ['y', 'yes']:
            break
        
        mutation = get_mutation_input(wt_handler)
        mutations.append(mutation)
    
    return mutations


def run_bioemu_generation(
    wt_sequence: SequenceHandler,
    mutation: Mutation,
    output_dir: Path,
    config: dict
) -> dict:
    """
    Run BioEmu to generate conformations.
    
    Args:
        wt_sequence: Wild-type sequence handler
        mutation: Mutation to study
        output_dir: Output directory
        config: Configuration dictionary
        
    Returns:
        Dictionary with generation results
    """
    print_header("BIOEMU CONFORMATION GENERATION")
    
    # Create mutation handler and get mutant sequence
    mutation_handler = MutationHandler(wt_sequence)
    mutant_handler = mutation_handler.get_mutant_handler(mutation)
    
    print_info(f"Wild-type: {wt_sequence.name}")
    print_info(f"Mutant: {mutant_handler.name}")
    print_info(f"Output directory: {output_dir}")
    print()
    
    # Initialize BioEmu
    bioemu_config = BioEmuConfig()
    if 'bioemu' in config:
        bioemu_config = BioEmuConfig(**config['bioemu'])
    
    print_info(f"Number of conformations: {bioemu_config.num_conformations}")
    print_info(f"Device: {bioemu_config.device}")
    print()
    
    runner = BioEmuRunner(bioemu_config)
    
    # Generate for both WT and mutant
    results = {}
    
    print(f"{Colors.BLUE}Generating wild-type conformations...{Colors.END}")
    wt_dir = output_dir / "wt"
    results['wt'] = runner.generate_conformations(
        str(wt_sequence),
        wt_dir,
        wt_sequence.name
    )
    print_success(f"WT generation complete: {wt_dir}")
    
    print()
    print(f"{Colors.BLUE}Generating mutant conformations...{Colors.END}")
    mutant_dir = output_dir / "mutant"
    results['mutant'] = runner.generate_conformations(
        str(mutant_handler),
        mutant_dir,
        mutant_handler.name
    )
    print_success(f"Mutant generation complete: {mutant_dir}")
    
    # Create standalone scripts for HPC execution
    print()
    print_info("Creating standalone BioEmu scripts for HPC execution...")
    
    wt_script = create_bioemu_script(
        str(wt_sequence), wt_dir, wt_sequence.name, bioemu_config
    )
    mutant_script = create_bioemu_script(
        str(mutant_handler), mutant_dir, mutant_handler.name, bioemu_config
    )
    
    print_success(f"WT script: {wt_script}")
    print_success(f"Mutant script: {mutant_script}")
    
    results['scripts'] = {
        'wt': str(wt_script),
        'mutant': str(mutant_script)
    }
    
    return results


def run_alignment(output_dir: Path, config: dict) -> dict:
    """
    Align generated conformations.
    
    Args:
        output_dir: Directory with generated conformations
        config: Configuration dictionary
        
    Returns:
        Alignment results
    """
    print_header("STRUCTURE ALIGNMENT")
    
    # Load alignment config
    align_config = AlignmentConfig()
    if 'alignment' in config:
        align_config = AlignmentConfig(**config['alignment'])
    
    print_info(f"Alignment method: {align_config.method}")
    print_info(f"Selection: {align_config.selection}")
    print_info(f"Reference: {align_config.reference}")
    print()
    
    results = {}
    
    for variant in ['wt', 'mutant']:
        variant_dir = output_dir / variant
        if not variant_dir.exists():
            print_warning(f"Skipping {variant}: directory not found")
            continue
        
        # Find PDB files
        pdb_files = list(variant_dir.glob("*_conf_*.pdb"))
        
        if not pdb_files:
            print_warning(f"No PDB files found in {variant_dir}")
            print_info("Note: BioEmu needs to be run first to generate conformations")
            continue
        
        print(f"{Colors.BLUE}Aligning {variant} conformations ({len(pdb_files)} structures)...{Colors.END}")
        
        aligned_dir = variant_dir / "aligned"
        aligned_dir.mkdir(exist_ok=True)
        
        # Use MDTraj if available
        try:
            from src.alignment import align_conformations_mdtraj
            results[variant] = align_conformations_mdtraj(
                pdb_files, aligned_dir, align_config.selection
            )
            print_success(f"{variant.upper()} alignment complete")
            print_info(f"RMSD: {results[variant]['rmsd_mean']:.2f} ± {results[variant]['rmsd_std']:.2f} Å")
        except Exception as e:
            print_warning(f"Could not align {variant}: {e}")
    
    return results


def display_analysis_options():
    """Display available analysis methods."""
    print_header("CONFORMATIONAL ANALYSIS OPTIONS")
    
    print("""
After generating and aligning conformations, you can perform the following analyses:

{0}1. STRUCTURAL METRICS{1}
   • RMSD (Root Mean Square Deviation): Measure structural similarity
   • RMSF (Root Mean Square Fluctuation): Per-residue flexibility
   • Radius of Gyration (Rg): Overall compactness

{0}2. CONTACT ANALYSIS{1}
   • Contact Maps: Residue-residue proximity
   • Native Contacts: Compare to reference structure
   • Differential Contact Maps: WT vs Mutant differences

{0}3. SECONDARY STRUCTURE{1}
   • DSSP Analysis: Helix/Sheet/Coil content over ensemble
   • Secondary Structure Propensity: Per-residue tendencies

{0}4. PRINCIPAL COMPONENT ANALYSIS (PCA){1}
   • Identify dominant motions
   • Project conformations onto principal components
   • Compare WT/Mutant conformational spaces

{0}5. FREE ENERGY LANDSCAPE{1}
   • 2D projections (PC1 vs PC2)
   • Identify metastable states
   • Compare energy basins between WT and Mutant

{0}6. CLUSTERING{1}
   • Group similar conformations
   • Identify representative structures
   • Compare cluster populations

{0}7. DYNAMIC NETWORK ANALYSIS{1}
   • Correlation networks
   • Identify allosteric pathways
   • Community detection

{0}8. BINDING SITE ANALYSIS (for binder design){1}
   • Identify conformations suitable for binder design
   • Analyze mutation-proximal regions
   • Assess binding pocket flexibility

{2}How to run analysis:{1}
   python scripts/analyze_conformations.py --input <output_dir>

{2}Required packages:{1}
   - MDTraj/MDAnalysis for trajectory analysis
   - ProDy for PCA and dynamics
   - NetworkX for network analysis
   - scikit-learn for clustering
""".format(Colors.BOLD, Colors.END, Colors.CYAN))


def save_run_summary(
    output_dir: Path,
    wt_sequence: SequenceHandler,
    mutation: Mutation,
    bioemu_results: dict,
    alignment_results: dict
) -> Path:
    """Save a summary of the pipeline run."""
    summary = {
        "timestamp": datetime.now().isoformat(),
        "wild_type": {
            "name": wt_sequence.name,
            "sequence": str(wt_sequence),
            "length": len(wt_sequence),
        },
        "mutation": str(mutation),
        "output_directory": str(output_dir),
        "bioemu": bioemu_results,
        "alignment": alignment_results,
    }
    
    summary_path = output_dir / "run_summary.json"
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    
    return summary_path


def main():
    """Main pipeline execution."""
    print_header("MUTATION CONFORMATION PIPELINE")
    print("This pipeline generates and analyzes protein conformations")
    print("for mutations using BioEmu deep learning model.\n")
    
    # Load configuration
    config = load_config()
    
    # Step 1: Get wild-type sequence
    wt_sequence = get_sequence_input()
    
    # Step 2: Get mutation
    mutation = get_mutation_input(wt_sequence)
    
    # Ask for additional mutations
    mutations = get_additional_mutations(wt_sequence, [mutation])
    
    if len(mutations) > 1:
        print()
        print_info(f"Total mutations to study: {len(mutations)}")
        for m in mutations:
            print(f"  • {m}")
    
    # Step 3: Set up output directory
    output_base = PROJECT_ROOT / "data" / "outputs"
    
    for mut in mutations:
        mutation_dir = create_mutation_directory(output_base, wt_sequence.name, mut)
        print()
        print_info(f"Output directory for {mut}: {mutation_dir}")
        
        # Step 4: Run BioEmu
        bioemu_results = run_bioemu_generation(wt_sequence, mut, mutation_dir, config)
        
        # Step 5: Alignment (if structures exist)
        alignment_results = run_alignment(mutation_dir, config)
        
        # Save summary
        summary_path = save_run_summary(
            mutation_dir, wt_sequence, mut, bioemu_results, alignment_results
        )
        print_success(f"Run summary saved: {summary_path}")
    
    # Step 6: Display analysis options
    display_analysis_options()
    
    print_header("PIPELINE COMPLETE")
    print_success("Conformations generated successfully!")
    print()
    print("Next steps:")
    print("  1. Run the generated BioEmu scripts on your HPC cluster")
    print("  2. Run the analysis script on the generated conformations")
    print("  3. Compare WT and mutant conformational ensembles")
    print("  4. Use insights for binder design")
    print()


if __name__ == "__main__":
    main()

import parmed as pmd
import argparse


def load_amber_files(prmtop_file, inpcrd_file, output_top, output_gro):
    """Load Amber topology and coordinate files."""
    amber = pmd.load_file(prmtop_file, inpcrd_file)
    amber.save(output_top, overwrite=True)
    amber.save(output_gro, overwrite=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Load Amber files")
    parser.add_argument("--prmtop", required=True, help="Path to .prmtop file")
    parser.add_argument("--inpcrd", required=True, help="Path to .inpcrd file")
    parser.add_argument("--output_top", required=True, help="Output .top file")
    parser.add_argument("--output_gro", required=True, help="Output .gro file")
    
    args = parser.parse_args()
    
    load_amber_files(args.prmtop, args.inpcrd, args.output_top, args.output_gro)

import os
import argparse
import sys
import subprocess

out_tc_groups = None
out_comm_grps = None

npt1 = """
define 			        = -DPOSRES
integrator              = md
dt                      = {dt} 
nsteps                  = {nsteps} 
nstxtcout               = {nstxtcout} 
nstvout                 = {nstvout} 
nstfout                 = {nstfout} 
nstcalcenergy           = {nstcalcenergy} 
nstenergy               = {nstenergy} 
nstlog                  = {nstlog} 
;
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 0.9
vdwtype                 = Cut-off
vdw-modifier            = None
DispCorr                = EnerPres
rvdw                    = 0.9
coulombtype             = PME
rcoulomb                = 0.9
;
tcoupl                  = berendsen
tc_grps                 = {tc_grps} 
tau_t                   = 1.0 1.0 
ref_t                   = 323.15 323.15 
;
pcoupl                  = berendsen
pcoupltype              = semiisotropic 
tau_p                   = 5.0
compressibility         = 4.5e-5  4.5e-5
ref_p                   = 1.0     1.0
refcoord_scaling        = com
;
constraints             = h-bonds
constraint_algorithm    = LINCS
continuation            = yes
;
nstcomm                 = 100
comm_mode               = linear
comm_grps               =  {comm_grps} 
"""


npt2 = """
define 			= -DPOSRES
integrator              = md
dt                      = {dt} 
nsteps                  = {nsteps} 
nstxtcout               = {nstxtcout} 
nstvout                 = {nstvout} 
nstfout                 = {nstfout} 
nstcalcenergy           = {nstcalcenergy} 
nstenergy               = {nstenergy} 
nstlog                  = {nstlog} 
;
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 0.9
vdwtype                 = Cut-off
vdw-modifier            = None
DispCorr                = EnerPres
rvdw                    = 0.9
coulombtype             = PME
rcoulomb                = 0.9
;
tcoupl                  = berendsen
tc_grps                 = {tc_grps}   
tau_t                   = 1.0 1.0 
ref_t                   = 323.15 323.15 
;
pcoupl                  = berendsen
pcoupltype              = semiisotropic 
tau_p                   = 5.0
compressibility         = 4.5e-5  4.5e-5
ref_p                   = 1.0     1.0
refcoord_scaling        = com
;
constraints             = h-bonds
constraint_algorithm    = LINCS
continuation            = yes
;
nstcomm                 = 100
comm_mode               = linear
comm_grps               = {comm_grps} 
"""

npt3 = """  
define 			= -DPOSRES
integrator              = md
dt                      = {dt} 
nsteps                  = {nsteps} 
nstxtcout               = {nstxtcout} 
nstvout                 = {nstvout} 
nstfout                 = {nstfout} 
nstcalcenergy           = {nstcalcenergy} 
nstenergy               = {nstenergy}
nstlog                  = {nstlog} 
;
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 0.9
vdwtype                 = Cut-off
vdw-modifier            = None
DispCorr                = EnerPres
rvdw                    = 0.9
coulombtype             = PME
rcoulomb                = 0.9
;
tcoupl                  = berendsen
tc_grps                 = {tc_grps} 
tau_t                   = 1.0 1.0 
ref_t                   = 323.15 323.15 
;
pcoupl                  = berendsen
pcoupltype              = semiisotropic 
tau_p                   = 5.0
compressibility         = 4.5e-5  4.5e-5
ref_p                   = 1.0     1.0
refcoord_scaling        = com
;
constraints             = h-bonds
constraint_algorithm    = LINCS
continuation            = yes
;
nstcomm                 = 100
comm_mode               = linear
comm_grps               = {comm_grps} 
"""


npt4 = """
define 			= -DPOSRES
integrator              = md
dt                      = {dt} 
nsteps                  = {nsteps} 
nstxtcout               = {nstxtcout} 
nstvout                 = {nstvout} 
nstfout                 = {nstfout} 
nstcalcenergy           = {nstcalcenergy} 
nstenergy               = {nstenergy} 
nstlog                  = {nstlog} 
;
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 0.9
vdwtype                 = Cut-off
vdw-modifier            = None
DispCorr                = EnerPres
rvdw                    = 0.9
coulombtype             = PME
rcoulomb                = 0.9
;
tcoupl                  = berendsen
tc_grps                 = {tc_grps} 
tau_t                   = 1.0 1.0
ref_t                   = 323.15 323.15 
;
pcoupl                  = berendsen
pcoupltype              = semiisotropic 
tau_p                   = 5.0
compressibility         = 4.5e-5  4.5e-5
ref_p                   = 1.0     1.0
refcoord_scaling        = com
;
constraints             = h-bonds
constraint_algorithm    = LINCS
continuation            = yes
;
nstcomm                 = 100
comm_mode               = linear
comm_grps               = {comm_grps} 
"""

nvt1 = """
define 			= -DPOSRES
integrator              = md
dt                      = {dt} 
nsteps                  = {nsteps} 
nstxtcout               = {nstxtcout} 
nstvout                 = {nstvout} 
nstfout                 = {nstfout} 
nstcalcenergy           = {nstcalcenergy} 
nstenergy               = {nstenergy} 
nstlog                  = {nstlog} 
;
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 0.9
vdwtype                 = Cut-off
vdw-modifier            = None
DispCorr                = EnerPres
rvdw                    = 0.9
coulombtype             = PME
rcoulomb                = 0.9
;
tcoupl                  = berendsen
tc_grps                 = {tc_grps} 
tau_t                   = 1.0 1.0
ref_t                   = 323.15 323.15
;
constraints             = h-bonds
constraint_algorithm    = LINCS
;
nstcomm                 = 100
comm_mode               = linear
comm_grps               = {comm_grps} 
;
gen-vel                 = yes
gen-temp                = 323.15
gen-seed                = -1
"""

nvt2 = """
define 			= -DPOSRES
integrator              = md
dt                      = {dt} 
nsteps                  = {nsteps} 
nstxtcout               = {nstxtcout} 
nstvout                 = {nstvout} 
nstfout                 = {nstfout} 
nstcalcenergy           = {nstcalcenergy} 
nstenergy               = {nstenergy} 
nstlog                  = {nstlog} 
;
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 0.9
vdwtype                 = Cut-off
vdw-modifier            = None
DispCorr                = EnerPres
rvdw                    = 0.9
coulombtype             = PME
rcoulomb                = 0.9
;
tcoupl                  = berendsen
tc_grps                 = {tc_grps} 
tau_t                   = 1.0 1.0
ref_t                   = 323.15 323.15
;
constraints             = h-bonds
constraint_algorithm    = LINCS
continuation            = yes
;
nstcomm                 = 100
comm_mode               = linear
comm_grps               = {comm_grps} 
"""

production = """    
; PRODUCTION MD RUN (10 ns, 5,000,000 steps at 2 fs)

integrator              = md
dt                      = {dt}          ; 2 fs
nsteps                  = {nsteps}        ; 10,000 ps = 10 ns

; Trajectory / output control
nstxout                 = {nstxout}              ; no full-precision coords to .trr
nstvout                 = {nstvout}
nstfout                 = {nstfout}
nstxout-compressed      = {nstxout_compressed}           ; coords every 10 ps to .xtc

nstcalcenergy           = {nstcalcenergy}          ; calc energies every 20 ps
nstenergy               = {nstenergy}          ; write energies every 20 ps
nstlog                  = {nstlog}          ; log every 20 ps

; Neighborsearching / cutoffs
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
DispCorr                = EnerPres
rvdw_switch             = 1.0
rvdw                    = 1.2
coulombtype             = PME
rcoulomb                = 1.2

; Temperature coupling (production)
tcoupl                  = Nose-Hoover
tc_grps                 = {tc_grps} 
tau_t                   = 1.0 1.0
ref_t                   = 323.15 323.15

; Pressure coupling (production)
pcoupl                  = Parrinello-Rahman
pcoupltype              = semiisotropic
tau_p                   = 5.0
compressibility         = 4.5e-5  4.5e-5
ref_p                   = 1.0     1.0

; Constraints
constraints             = h-bonds
constraint_algorithm    = LINCS
continuation            = yes

; Center-of-mass motion removal
nstcomm                 = 100
comm_mode               = linear
comm_grps               =  {comm_grps} 
"""

# Default parameter values for each template
DEFAULTS = {
    'npt1': {
        'dt': 0.001,
        'nsteps': 125000,
        'nstxtcout': 5000,
        'nstvout': 5000,
        'nstfout': 5000,
        'nstcalcenergy': 100,
        'nstenergy': 1000,
        'nstlog': 1000,
        'tc_grps': out_tc_groups,
        'comm_grps': out_comm_grps
    },
    'npt2': {
        'dt': 0.002,
        'nsteps': 250000,
        'nstxtcout': 5000,
        'nstvout': 5000,
        'nstfout': 5000,
        'nstcalcenergy': 100,
        'nstenergy': 1000,
        'nstlog': 1000,
        'tc_grps': out_tc_groups,
        'comm_grps': out_comm_grps
    },
    'npt3': {
        'dt': 0.002,
        'nsteps': 250000,
        'nstxtcout': 5000,
        'nstvout': 5000,
        'nstfout': 5000,
        'nstcalcenergy': 100,
        'nstenergy': 1000,
        'nstlog': 1000,
        'tc_grps': out_tc_groups,
        'comm_grps': out_comm_grps
    },
    'npt4': {
        'dt': 0.002,
        'nsteps': 250000,
        'nstxtcout': 5000,
        'nstvout': 5000,
        'nstfout': 5000,
        'nstcalcenergy': 100,
        'nstenergy': 1000,
        'nstlog': 1000,
        'tc_grps': out_tc_groups,
        'comm_grps': out_comm_grps
    },
    'nvt1': {
        'dt': 0.001,
        'nsteps': 125000,
        'nstxtcout': 5000,
        'nstvout': 5000,
        'nstfout': 5000,
        'nstcalcenergy': 100,
        'nstenergy': 1000,
        'nstlog': 1000,
        'tc_grps': out_tc_groups,
        'comm_grps': out_comm_grps
    },
    'nvt2': {
        'dt': 0.001,
        'nsteps': 125000,
        'nstxtcout': 5000,
        'nstvout': 5000,
        'nstfout': 5000,
        'nstcalcenergy': 100,
        'nstenergy': 1000,
        'nstlog': 1000,
        'tc_grps': out_tc_groups,
        'comm_grps': out_comm_grps
    },
    'production': {
        'dt': 0.002,
        'nsteps': 5000000,
        'nstxout': 0,
        'nstvout': 0,
        'nstfout': 0,
        'nstxout_compressed': 10000,
        'nstcalcenergy': 10000,
        'nstenergy': 10000,
        'nstlog': 10000,
        'tc_grps': out_tc_groups,
        'comm_grps': out_comm_grps
    }
}

TEMPLATE_MAP = {
    'npt1': npt1,
    'npt2': npt2,
    'npt3': npt3,
    'npt4': npt4,
    'nvt1': nvt1,
    'nvt2': nvt2,
    'production': production
}


def create_parser():
    """Create argument parser with subparsers for each template."""
    parser = argparse.ArgumentParser(
        description='Generate GROMACS MDP files from templates'
    )
    
    subparsers = parser.add_subparsers(dest='template', help='Template to use')
    
    # Define common parameters
    common_params = {
        '--dt': {'type': float, 'help': 'Time step (ps)'},
        '--nsteps': {'type': int, 'help': 'Number of steps'},
        '--nstxtcout': {'type': int, 'help': 'XTC output frequency'},
        '--nstvout': {'type': int, 'help': 'Velocity output frequency'},
        '--nstfout': {'type': int, 'help': 'Force output frequency'},
        '--nstcalcenergy': {'type': int, 'help': 'Energy calculation frequency'},
        '--nstenergy': {'type': int, 'help': 'Energy output frequency'},
        '--nstlog': {'type': int, 'help': 'Log output frequency'},
        '--tc_grps': {'type': str, 'help': 'Temperature coupling groups'},
        '--comm_grps': {'type': str, 'help': 'COM motion removal groups'}
    }
    
    production_params = {
        '--dt': {'type': float, 'help': 'Time step (ps)'},
        '--nsteps': {'type': int, 'help': 'Number of steps'},
        '--nstxout': {'type': int, 'help': 'Full precision output frequency'},
        '--nstvout': {'type': int, 'help': 'Velocity output frequency'},
        '--nstfout': {'type': int, 'help': 'Force output frequency'},
        '--nstxout_compressed': {'type': int, 'help': 'Compressed trajectory output frequency'},
        '--nstcalcenergy': {'type': int, 'help': 'Energy calculation frequency'},
        '--nstenergy': {'type': int, 'help': 'Energy output frequency'},
        '--nstlog': {'type': int, 'help': 'Log output frequency'},
        '--tc_grps': {'type': str, 'help': 'Temperature coupling groups'},
        '--comm_grps': {'type': str, 'help': 'COM motion removal groups'}
    }
    
    # Define custom parameters (all options from generate_mdp.py)
    custom_params = {
        # Run control
        '--integrator': {'type': str, 'help': 'Integrator type (e.g., md, steep)'},
        '--dt': {'type': str, 'help': 'Time step (ps)'},
        '--nsteps': {'type': str, 'help': 'Number of steps'},
        '--tinit': {'type': str, 'help': 'Initial simulation time'},
        '--init_step': {'type': str, 'help': 'Initial step number'},
        '--continuation': {'type': str, 'help': 'Continuation of simulation (yes/no)'},
        # Neighbor searching
        '--cutoff_scheme': {'type': str, 'help': 'Cutoff scheme (Verlet/Group)'},
        '--nstlist': {'type': str, 'help': 'Neighbor list update frequency'},
        '--pbc': {'type': str, 'help': 'Periodic boundary conditions (xyz/xy/scxy)'},
        '--verlet_buffer_tolerance': {'type': str, 'help': 'Verlet buffer tolerance'},
        '--rlist': {'type': str, 'help': 'Neighbor list cutoff'},
        '--rcoulomb': {'type': str, 'help': 'Coulomb cutoff'},
        '--rvdw': {'type': str, 'help': 'VDW cutoff'},
        # Electrostatics and vdw
        '--coulombtype': {'type': str, 'help': 'Coulomb interaction type (PME/Reaction-Field/Cut-off)'},
        '--coulomb_modifier': {'type': str, 'help': 'Coulomb modifier (Potential-shift)'},
        '--vdwtype': {'type': str, 'help': 'VDW interaction type (Cut-off/Shift/Switch)'},
        '--vdw_modifier': {'type': str, 'help': 'VDW modifier (Potential-shift/Force-switch)'},
        '--dispcorr': {'type': str, 'help': 'Dispersion correction (EnerPres/Ener/No)'},
        '--fourierspacing': {'type': str, 'help': 'Fourier spacing'},
        '--pme_order': {'type': str, 'help': 'PME interpolation order'},
        '--ewald_rtol': {'type': str, 'help': 'Ewald relative tolerance'},
        # Temperature coupling
        '--tcoupl': {'type': str, 'help': 'Temperature coupling algorithm (v-rescale/Nose-Hoover/berendsen)'},
        '--tc_grps': {'type': str, 'help': 'Temperature coupling groups'},
        '--tau_t': {'type': str, 'help': 'Temperature coupling time constants'},
        '--ref_t': {'type': str, 'help': 'Reference temperatures'},
        '--nsttcouple': {'type': str, 'help': 'Temperature coupling frequency'},
        # Velocity generation
        '--gen_vel': {'type': str, 'help': 'Generate velocities (yes/no)'},
        '--gen_temp': {'type': str, 'help': 'Temperature for velocity generation'},
        '--gen_seed': {'type': str, 'help': 'Seed for velocity generation'},
        # Pressure coupling
        '--pcoupl': {'type': str, 'help': 'Pressure coupling (Parrinello-Rahman/C-rescale/berendsen/no)'},
        '--pcoupltype': {'type': str, 'help': 'Pressure coupling type (isotropic/semiisotropic/anisotropic)'},
        '--tau_p': {'type': str, 'help': 'Pressure coupling time constant'},
        '--ref_p': {'type': str, 'help': 'Reference pressure'},
        '--compressibility': {'type': str, 'help': 'Compressibility'},
        '--nstpcouple': {'type': str, 'help': 'Pressure coupling frequency'},
        # Constraints
        '--constraints': {'type': str, 'help': 'Constraint algorithm (h-bonds/all-bonds/h-angles/all-angles)'},
        '--constraint_algorithm': {'type': str, 'help': 'Constraint algorithm (LINCS/SHAKE)'},
        '--lincs_order': {'type': str, 'help': 'LINCS order'},
        '--lincs_iter': {'type': str, 'help': 'LINCS iterations'},
        '--shake_tol': {'type': str, 'help': 'SHAKE tolerance'},
        # Output control
        '--nstxout': {'type': str, 'help': 'Full precision output frequency'},
        '--nstvout': {'type': str, 'help': 'Velocity output frequency'},
        '--nstfout': {'type': str, 'help': 'Force output frequency'},
        '--nstlog': {'type': str, 'help': 'Log output frequency'},
        '--nstenergy': {'type': str, 'help': 'Energy output frequency'},
        '--nstxout_compressed': {'type': str, 'help': 'Compressed trajectory output frequency'},
        '--compressed_x_precision': {'type': str, 'help': 'Compressed coordinate precision'},
        # Pulling
        '--pull': {'type': str, 'help': 'Pulling (yes/no)'},
        '--pull_ngroups': {'type': str, 'help': 'Number of pull groups'},
        '--pull_ncoords': {'type': str, 'help': 'Number of pull coordinates'},
        '--pull_groups': {'type': str, 'help': 'Pull groups'},
        '--pull_coord1_type': {'type': str, 'help': 'Pull coordinate type'},
        '--pull_coord1_geometry': {'type': str, 'help': 'Pull coordinate geometry'},
        '--pull_coord1_dim': {'type': str, 'help': 'Pull coordinate dimensions'},
        '--pull_coord1_vec': {'type': str, 'help': 'Pull coordinate vector'},
        '--pull_coord1_init': {'type': str, 'help': 'Pull coordinate initial value'},
        '--pull_coord1_rate': {'type': str, 'help': 'Pull coordinate rate'},
        '--pull_coord1_k': {'type': str, 'help': 'Pull coordinate force constant'},
        # Extra parameters
        '--extra_mdp': {'type': str, 'help': 'Extra MDP lines'},
        '--output': {'type': str, 'help': 'Output file name'}
    }
    
    # Create subparsers for each template
    for template_name in ['npt1', 'npt2', 'npt3', 'npt4', 'nvt1', 'nvt2']:
        sp = subparsers.add_parser(template_name, help=f'{template_name} template')
        for param, kwargs in common_params.items():
            sp.add_argument(param, **kwargs)
    
    # Production template with different parameters
    sp = subparsers.add_parser('production', help='Production template')
    for param, kwargs in production_params.items():
        sp.add_argument(param, **kwargs)
    
    # Index file builder
    sp = subparsers.add_parser('index', help='Generate GROMACS index file')
    sp.add_argument('--config', type=str, required=True, help='Configuration file with lipid, ligand, and ion components')
    sp.add_argument('--index_log', type=str, required=True, help='index_log.txt file with index group definitions')
    sp.add_argument('--gro-file', type=str, required=True, help='Input GRO structure file')
    sp.add_argument('--output', type=str, required=True, help='Output index file name')
    
    # Custom template with all parameters from generate_mdp.py
    sp = subparsers.add_parser('custom', help='Custom MDP generation with full parameter control')
    for param, kwargs in custom_params.items():
        sp.add_argument(param, **kwargs)
    
    return parser


def merge_args_with_defaults(args, defaults):
    """Merge command-line arguments with default values."""
    params = defaults.copy()
    
    # Get all arguments from the namespace
    args_dict = vars(args)
    
    # Update with provided arguments (skip None values)
    for key, value in args_dict.items():
        if value is not None and key != 'template':
            params[key] = value
    
    return params


def format_template(template_str, params):
    """Format template with given parameters."""
    return template_str.format(**params)


def group_args_by_template(argv):
    """Group arguments by template name for multi-template processing."""
    template_names = {'npt1', 'npt2', 'npt3', 'npt4', 'nvt1', 'nvt2', 'production', 'index', 'custom'}
    groups = []
    current_template = None
    current_args = []
    
    for arg in argv:
        if arg in template_names:
            if current_template is not None:
                groups.append((current_template, current_args))
            current_template = arg
            current_args = []
        else:
            if current_template is not None:
                current_args.append(arg)
    
    if current_template is not None:
        groups.append((current_template, current_args))
    
    return groups


def process_single_template():
    """Process a single template with standard argparse."""
    parser = create_parser()
    args = parser.parse_args()
    
    if not args.template:
        parser.print_help()
        sys.exit(1)
    
    # Handle custom template differently
    if args.template == 'custom':
        try:
            # For custom template, all parameters come directly from command line
            params = {}
            args_dict = vars(args)
            
            for key, value in args_dict.items():
                if value is not None and key != 'template':
                    params[key] = value
            
            # Build MDP file content from parameters
            mdp_lines = []
            for key, value in params.items():
                if key != 'output' and key != 'extra_mdp' and value:
                    mdp_lines.append(f"{key:<22} = {value}")
            
            # Add extra MDP lines if provided
            if 'extra_mdp' in params and params['extra_mdp']:
                mdp_lines.append("")
                mdp_lines.append("; ---- Extra .mdp lines ----")
                for line in params['extra_mdp'].strip().split('\n'):
                    if line.strip():
                        mdp_lines.append(line)
            
            mdp_content = "\n".join(mdp_lines)
            
            # Write to output file or print
            output_file = params.get('output')
            if output_file:
                with open(output_file, 'w') as f:
                    f.write(mdp_content)
                print(f"Generated: {output_file}")
            else:
                print(mdp_content)
            
        except Exception as e:
            print(f"Error processing custom template: {e}", file=sys.stderr)
            sys.exit(1)
        return
    
    # Get defaults and merge with provided arguments
    defaults = DEFAULTS[args.template]
    params = merge_args_with_defaults(args, defaults)
    
    # Format the template
    try:
        template_str = TEMPLATE_MAP[args.template]
        formatted_mdp = format_template(template_str, params)
        print(formatted_mdp)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

def process_multiple_templates(argv):
    """Process multiple templates from command line arguments."""
    parser = create_parser()
    
    # Group arguments by template
    groups = group_args_by_template(argv)
    
    if not groups:
        print("Error: No templates specified", file=sys.stderr)
        sys.exit(1)
    
    # Check for custom template mutual exclusivity
    template_names = [g[0] for g in groups]
    if 'custom' in template_names:
        if len(template_names) > 1:
            print("Error: 'custom' template cannot be used with other templates (npt1, npt2, npt3, npt4, nvt1, nvt2, production, index)", file=sys.stderr)
            sys.exit(1)
        # Process custom template
        custom_template_args = groups[0][1]
        try:
            args = parser.parse_args(['custom'] + custom_template_args)
            
            # For custom template, all parameters come directly from command line
            params = {}
            args_dict = vars(args)
            
            for key, value in args_dict.items():
                if value is not None and key != 'template':
                    params[key] = value
            
            # Build MDP file content from parameters
            mdp_lines = []
            for key, value in params.items():
                if key != 'output' and key != 'extra_mdp' and value:
                    mdp_lines.append(f"{key:<22} = {value}")
            
            # Add extra MDP lines if provided
            if 'extra_mdp' in params and params['extra_mdp']:
                mdp_lines.append("")
                mdp_lines.append("; ---- Extra .mdp lines ----")
                for line in params['extra_mdp'].strip().split('\n'):
                    if line.strip():
                        mdp_lines.append(line)
            
            mdp_content = "\n".join(mdp_lines)
            
            # Write to output file
            output_file = params.get('output')
            if not output_file:
                print("Error: --output is required for custom template", file=sys.stderr)
                sys.exit(1)
            
            with open(output_file, 'w') as f:
                f.write(mdp_content)
            print(f"Generated: {output_file}")
            
        except SystemExit as e:
            if e.code != 0:
                sys.exit(1)
        except Exception as e:
            print(f"Error processing custom template: {e}", file=sys.stderr)
            sys.exit(1)
        return
    
    # Process standard templates (npt1, npt2, npt3, npt4, nvt1, nvt2, production, index)
    all_outputs = {}
    
    # Process each template group
    for template_name, template_args in groups:
        try:
            # Handle index template separately
            if template_name == 'index':
                args = parser.parse_args([template_name] + template_args)
                tc_groups, comm_grps = GenerateCommandForIndexFile(args)
                continue
            
            # Parse arguments for this template
            args = parser.parse_args([template_name] + template_args)
            
            if not args.template:
                continue
            
            # Get defaults and merge with provided arguments

            DEFAULTS['npt1']['tc_grps'] = tc_groups
            DEFAULTS['npt1']['comm_grps'] = comm_grps
            DEFAULTS['npt2']['tc_grps'] = tc_groups
            DEFAULTS['npt2']['comm_grps'] = comm_grps
            DEFAULTS['npt3']['tc_grps'] = tc_groups
            DEFAULTS['npt3']['comm_grps'] = comm_grps
            DEFAULTS['npt4']['tc_grps'] = tc_groups
            DEFAULTS['npt4']['comm_grps'] = comm_grps
            DEFAULTS['nvt1']['tc_grps'] = tc_groups
            DEFAULTS['nvt1']['comm_grps'] = comm_grps
            DEFAULTS['nvt2']['tc_grps'] = tc_groups
            DEFAULTS['nvt2']['comm_grps'] = comm_grps
            DEFAULTS['production']['tc_grps'] = tc_groups
            DEFAULTS['production']['comm_grps'] = comm_grps 

            defaults = DEFAULTS[args.template]

            params = merge_args_with_defaults(args, defaults)
            
            # Format the template
            template_str = TEMPLATE_MAP[args.template]
            formatted_mdp = format_template(template_str, params)
            all_outputs[template_name] = formatted_mdp
            
        except SystemExit:
            continue
        except ValueError as e:
            print(f"Error processing {template_name}: {e}", file=sys.stderr)
            sys.exit(1)
    
    for template_name, output in all_outputs.items():
        filename = f"{template_name}.mdp"
        try:
            with open(filename, 'w') as f:
                f.write(output)
            print(f"Generated: {filename}")
        except IOError as e:
            print(f"Error writing {filename}: {e}", file=sys.stderr)
            sys.exit(1)


def GenerateCommandForIndexFile(args):

    monovalent_ions = [
    "Na+", "NA+",
    "K+",  "K+",
    "Li+", "LI+",
    "Rb+", "RB+",
    "Cs+", "CS+",
    "Cl-", "CL-",
    "F-",  "F-",
    "Br-", "BR-",
    "I-",  "I-",
    ]


    tc_grps = []
    comm_grps = []

    index_dict = {}
    extracted_lines = []

    lipid_frgs = []
    ions = []
    ligand = []

    with open(args.config, 'r') as f:
        conf_lines = f.readlines()

    with open(args.index_log, 'r') as f:
        index_log_lines = f.readlines()

    for l in conf_lines:
        if "lipid"  in l:
            lipid_frgs = l.strip('\n').split(":")[1].split(',')
        elif "ligand" in l:
            ligand = l.strip('\n').split(":")[1].split(',')

    #ions extracted from conf files sometimes have different naming and therefore causing the index name missmatch. Hence, extracting ions from index log file and matching with monovalent ion list to get the correct ion names for index group creation.
    # elif "ions" in l:
    #     ions = l.strip('\n').split(":")[1].split(',') 

    for i in index_log_lines:
        if len(i.split())  == 5 and i.split()[2] == ":" and i.split()[4] == 'atoms':
            extracted_lines.append(i.split())
            index_dict[i.split()[1]] = i.split()[0]
            # print(i.split()[0], "-->", i.split()[1].strip())
            if i.split()[1].strip() in monovalent_ions:
                ions.append(i.split()[1].strip())
                
    lipid_grp_commands = []
    lipid_index = []
    inons_index = []

    end_ndx  = int(extracted_lines[-1][0])

    # print(extracted_lines)
    # lipid_frgs = ['PA', 'PC'] 
    # ions = ['NA+', 'CL-']
    # ligand = []

    inons_index = []
    for f in ions:
        inons_index.append(index_dict[f])

    lipid_index = []
    for f in lipid_frgs:
        lipid_index.append(index_dict[f])

    lipid_group_commands = []
    ligand_group_commands = []
    ions_group_commands = []

    lipid_protein_index = None
    new_index_grops = {}

    if len(lipid_frgs) == 0:
        pass 
    elif len(lipid_frgs) == 1:
        lipid_group_commands.append("echo '"+index_dict["Protein"]+" | "+index_dict[lipid_frgs[0]]+"'")
        end_ndx = end_ndx+1 
        lipid_protein_index = end_ndx
        new_index_grops[lipid_protein_index] = f"Protein_{lipid_frgs[0]}"

    elif len(lipid_frgs) > 1:
        lipid_group_commands.append(f"echo '{' | '.join(lipid_index)}'")
    
        end_ndx = end_ndx+1 
        lipid_group_commands.append("echo '"+index_dict["Protein"]+" | "+str(end_ndx) +"'")
        
        end_ndx = end_ndx + 1
        lipid_protein_index = end_ndx
        new_index_grops[lipid_protein_index] = f"Protein_{'_'.join(lipid_frgs)}"

        tc_grps.append(f"{'_'.join(lipid_frgs)}")

    if len(ligand) == 1:
        ligand_group_commands.append("echo '"+index_dict["Protein"]+ "|"+ index_dict[ligand[0]]+"'")
        end_ndx = end_ndx+1 
        new_index_grops[end_ndx] = f"Protein_{ligand[0]}"

        tc_grps.append(f"Protein_{ligand[0]}")
        
        if lipid_protein_index:
            ligand_group_commands.append("echo '"+str(lipid_protein_index)+ "|"+ index_dict[ligand[0]]+"'")
            end_ndx = end_ndx+1 
            new_index_grops[end_ndx] = f"{new_index_grops[lipid_protein_index]}_{ligand[0]}"

            comm_grps.append(f"{new_index_grops[lipid_protein_index]}_{ligand[0]}")
    else:
        tc_grps.append("Protein")
        comm_grps.append(f"{new_index_grops[lipid_protein_index]}")
        
    
    if len(ions) == 0:
        pass
    elif len(ions) == 1:
        ions_group_commands.append("echo '"+index_dict["Water"]+" | "+index_dict[ions[0]]+"'")
        end_ndx = end_ndx+1 
        new_index_grops[end_ndx] = f"Water_{ions[0]}"
        tc_grps.append(f"Water_{ions[0]}")
        comm_grps.append(f"Water_{ions[0]}")
        
    elif len(ions) > 1:
        ions_group_commands.append(f"echo '{' | '.join(inons_index)}'")
        end_ndx = end_ndx+1 
        ions_group_commands.append("echo '"+index_dict["Water"]+" | "+str(end_ndx) +"'")
        end_ndx = end_ndx + 1
        new_index_grops[end_ndx] = f"Water_{'_'.join(ions)}"
        tc_grps.append(f"Water_{'_'.join(ions)}")
        comm_grps.append(f"Water_{'_'.join(ions)}")

    command_list = []

    if len(lipid_group_commands) > 0:
        command_list.append(' ; '.join(lipid_group_commands))
        
    if len(ligand_group_commands) > 0:
        command_list.append(' ; '.join(ligand_group_commands))

    if len(ions_group_commands) > 0:
        command_list.append(' ; '.join(ions_group_commands))

 
    #DEBUGGING PRINTS       
    # print (f"({' ; '.join(command_list)} ; echo 'q') | gmx make_ndx -f system.gro -o ndx.ndx")
    # print(" ".join(tc_grps))
    # print(" ".join(comm_grps))
        
    cmd_string = f"({' ; '.join(command_list)} ; echo 'q') | gmx make_ndx -f system.gro -o ndx.ndx"

    # Execute the command
    try:
        result = subprocess.run(cmd_string, shell=True, check=True)
        print(f"\nIndex file successfully created: {args.output}")
    except subprocess.CalledProcessError as e:
        print(f"Error executing gmx make_ndx: {e}", file=sys.stderr)
        sys.exit(1)
    return " ".join(comm_grps), " ".join(comm_grps)
    
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python md_template.py [template] [options] [[template] [options] ...]")
        print("       python md_template.py [template] [options] --write  (write to file)")
        print("       python md_template.py index --config CONFIG --gro-file GRO --output OUTPUT")
        print("       python md_template.py npt1 [options] npt2 [options] index --config CONFIG --gro-file GRO --output OUTPUT")
        print("       python md_template.py custom [options] --output OUTPUT  (custom parameters)")
        print("Run 'python md_template.py npt1 -h' for help on a specific template")
        print("Note: 'custom' template cannot be used with other templates")
        sys.exit(1)
    
    # Remove write flag from argv before processing
    argv_clean = [arg for arg in sys.argv[1:] if arg not in ['-w', '--write']]
    
    # Check if multiple templates are requested (including index and custom)
    template_names = {'npt1', 'npt2', 'npt3', 'npt4', 'nvt1', 'nvt2', 'production', 'index', 'custom'}
    template_count = sum(1 for arg in argv_clean if arg in template_names)
    
    if template_count > 1:
        # Multiple templates mode
        process_multiple_templates(argv_clean)
    else:
        # Single template mode
        process_single_template()

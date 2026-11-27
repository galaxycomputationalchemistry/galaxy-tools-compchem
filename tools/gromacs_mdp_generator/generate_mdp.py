#!/usr/bin/env python3
import argparse
import re
from collections import OrderedDict

PRESETS = {
    "em": OrderedDict([
        ("integrator", "steep"),
        ("emtol", "1000.0"),
        ("emstep", "0.01"),
        ("nsteps", "50000"),
        ("cutoff-scheme", "Verlet"),
        ("nstlist", "10"),
        ("pbc", "xyz"),
        ("coulombtype", "PME"),
        ("rcoulomb", "1.0"),
        ("vdwtype", "Cut-off"),
        ("rvdw", "1.0"),
        ("constraints", "h-bonds"),
        ("constraint-algorithm", "LINCS"),
        ("lincs-order", "4"),
        ("lincs-iter", "1"),
        ("nstlog", "100"),
        ("nstenergy", "100"),
    ]),
    "nvt": OrderedDict([
        ("integrator", "md"),
        ("dt", "0.002"),
        ("nsteps", "50000"),
        ("continuation", "no"),
        ("cutoff-scheme", "Verlet"),
        ("nstlist", "20"),
        ("pbc", "xyz"),
        ("coulombtype", "PME"),
        ("coulomb-modifier", "Potential-shift"),
        ("rcoulomb", "1.0"),
        ("vdwtype", "Cut-off"),
        ("vdw-modifier", "Potential-shift"),
        ("rvdw", "1.0"),
        ("tcoupl", "v-rescale"),
        ("tc-grps", "System"),
        ("tau-t", "0.1"),
        ("ref-t", "300"),
        ("nsttcouple", "-1"),
        ("pcoupl", "no"),
        ("constraints", "h-bonds"),
        ("constraint-algorithm", "LINCS"),
        ("lincs-order", "4"),
        ("lincs-iter", "1"),
        ("gen-vel", "yes"),
        ("gen-temp", "300"),
        ("gen-seed", "-1"),
        ("nstlog", "1000"),
        ("nstenergy", "1000"),
        ("nstxout-compressed", "1000"),
        ("compressed-x-precision", "1000"),
    ]),
    "npt": OrderedDict([
        ("integrator", "md"),
        ("dt", "0.002"),
        ("nsteps", "50000"),
        ("continuation", "yes"),
        ("cutoff-scheme", "Verlet"),
        ("nstlist", "20"),
        ("pbc", "xyz"),
        ("coulombtype", "PME"),
        ("coulomb-modifier", "Potential-shift"),
        ("rcoulomb", "1.0"),
        ("vdwtype", "Cut-off"),
        ("vdw-modifier", "Potential-shift"),
        ("rvdw", "1.0"),
        ("tcoupl", "v-rescale"),
        ("tc-grps", "System"),
        ("tau-t", "0.1"),
        ("ref-t", "300"),
        ("nsttcouple", "-1"),
        ("pcoupl", "C-rescale"),
        ("pcoupltype", "isotropic"),
        ("tau-p", "5.0"),
        ("ref-p", "1.0"),
        ("compressibility", "4.5e-5"),
        ("nstpcouple", "-1"),
        ("constraints", "h-bonds"),
        ("constraint-algorithm", "LINCS"),
        ("lincs-order", "4"),
        ("lincs-iter", "1"),
        ("nstlog", "1000"),
        ("nstenergy", "1000"),
        ("nstxout-compressed", "1000"),
        ("compressed-x-precision", "1000"),
    ]),
    "md": OrderedDict([
        ("integrator", "md"),
        ("dt", "0.002"),
        ("nsteps", "250000"),
        ("continuation", "yes"),
        ("cutoff-scheme", "Verlet"),
        ("nstlist", "20"),
        ("pbc", "xyz"),
        ("coulombtype", "PME"),
        ("coulomb-modifier", "Potential-shift"),
        ("rcoulomb", "1.0"),
        ("vdwtype", "Cut-off"),
        ("vdw-modifier", "Potential-shift"),
        ("rvdw", "1.0"),
        ("tcoupl", "v-rescale"),
        ("tc-grps", "System"),
        ("tau-t", "0.5"),
        ("ref-t", "300"),
        ("nsttcouple", "-1"),
        ("pcoupl", "Parrinello-Rahman"),
        ("pcoupltype", "isotropic"),
        ("tau-p", "5.0"),
        ("ref-p", "1.0"),
        ("compressibility", "4.5e-5"),
        ("nstpcouple", "-1"),
        ("constraints", "h-bonds"),
        ("constraint-algorithm", "LINCS"),
        ("lincs-order", "4"),
        ("lincs-iter", "1"),
        ("nstlog", "1000"),
        ("nstenergy", "1000"),
        ("nstxout-compressed", "1000"),
        ("compressed-x-precision", "1000"),
    ]),
    "blank": OrderedDict([]),
}

AUTO = "__AUTO__"

def parse_extra_lines(extra: str):
    """
    Parse user-provided mdp lines. Keeps raw lines, but also extracts key/value.
    Supports:
      key = value
      key=value
      ; comments
      blank lines
    """
    lines = []
    kv = OrderedDict()
    if not extra:
        return lines, kv
    for raw in extra.splitlines():
        line = raw.rstrip("\n")
        lines.append(line)
        s = line.strip()
        if not s or s.startswith(";"):
            continue
        m = re.match(r"^([A-Za-z0-9_.-]+)\s*=\s*(.*)$", s)
        if m:
            k = m.group(1).strip()
            v = m.group(2).strip()
            kv[k] = v
    return lines, kv

def set_if_not_auto(params: OrderedDict, key: str, value: str):
    if value is None:
        return
    if value == AUTO:
        return
    # Allow blank to mean "don't set"
    if isinstance(value, str) and value.strip() == "":
        return
    params[key] = str(value)

def main():
    ap = argparse.ArgumentParser(description="Generate a GROMACS .mdp file (professional preset + overrides).")

    ap.add_argument("--preset", choices=list(PRESETS.keys()), required=True)

    # Run control
    ap.add_argument("--integrator", default="")
    ap.add_argument("--dt", default="")
    ap.add_argument("--nsteps", default="")
    ap.add_argument("--tinit", default="")
    ap.add_argument("--init_step", default="")
    ap.add_argument("--continuation", default="")

    # Nonbonded/neighbor
    ap.add_argument("--cutoff_scheme", default="")
    ap.add_argument("--nstlist", default="")
    ap.add_argument("--pbc", default="")
    ap.add_argument("--verlet_buffer_tolerance", default="")
    ap.add_argument("--rlist", default="")
    ap.add_argument("--rcoulomb", default="")
    ap.add_argument("--rvdw", default="")

    # Electrostatics + vdw
    ap.add_argument("--coulombtype", default="")
    ap.add_argument("--coulomb_modifier", default="")
    ap.add_argument("--vdwtype", default="")
    ap.add_argument("--vdw_modifier", default="")
    ap.add_argument("--dispcorr", default="")
    ap.add_argument("--fourierspacing", default="")
    ap.add_argument("--pme_order", default="")
    ap.add_argument("--ewald_rtol", default="")

    # Temperature / velocity generation
    ap.add_argument("--tcoupl", default="")
    ap.add_argument("--tc_grps", default="")
    ap.add_argument("--tau_t", default="")
    ap.add_argument("--ref_t", default="")
    ap.add_argument("--nsttcouple", default="")
    ap.add_argument("--gen_vel", default="")
    ap.add_argument("--gen_temp", default="")
    ap.add_argument("--gen_seed", default="")

    # Pressure coupling
    ap.add_argument("--pcoupl", default="")
    ap.add_argument("--pcoupltype", default="")
    ap.add_argument("--tau_p", default="")
    ap.add_argument("--ref_p", default="")
    ap.add_argument("--compressibility", default="")
    ap.add_argument("--nstpcouple", default="")

    # Constraints
    ap.add_argument("--constraints", default="")
    ap.add_argument("--constraint_algorithm", default="")
    ap.add_argument("--lincs_order", default="")
    ap.add_argument("--lincs_iter", default="")
    ap.add_argument("--shake_tol", default="")

    # Output
    ap.add_argument("--nstxout", default="")
    ap.add_argument("--nstvout", default="")
    ap.add_argument("--nstfout", default="")
    ap.add_argument("--nstlog", default="")
    ap.add_argument("--nstenergy", default="")
    ap.add_argument("--nstxout_compressed", default="")
    ap.add_argument("--compressed_x_precision", default="")

    # Pull (basic)
    ap.add_argument("--pull", default=AUTO)
    ap.add_argument("--pull_ngroups", default="")
    ap.add_argument("--pull_ncoords", default="")
    ap.add_argument("--pull_groups", default="")
    ap.add_argument("--pull_coord1_type", default="")
    ap.add_argument("--pull_coord1_geometry", default="")
    ap.add_argument("--pull_coord1_dim", default="")
    ap.add_argument("--pull_coord1_vec", default="")
    ap.add_argument("--pull_coord1_init", default="")
    ap.add_argument("--pull_coord1_rate", default="")
    ap.add_argument("--pull_coord1_k", default="")

    # Expert overrides
    ap.add_argument("--extra_mdp", default="")

    ap.add_argument("--output", required=True)

    args = ap.parse_args()

    params = OrderedDict(PRESETS.get(args.preset, OrderedDict()))

    # Map UI args to mdp keys
    set_if_not_auto(params, "integrator", args.integrator)
    set_if_not_auto(params, "continuation", args.continuation)
    set_if_not_auto(params, "cutoff-scheme", args.cutoff_scheme)
    set_if_not_auto(params, "pbc", args.pbc)
    set_if_not_auto(params, "coulombtype", args.coulombtype)
    set_if_not_auto(params, "coulomb-modifier", args.coulomb_modifier)
    set_if_not_auto(params, "vdwtype", args.vdwtype)
    set_if_not_auto(params, "vdw-modifier", args.vdw_modifier)
    set_if_not_auto(params, "DispCorr", args.dispcorr)
    set_if_not_auto(params, "tcoupl", args.tcoupl)
    set_if_not_auto(params, "pcoupl", args.pcoupl)
    set_if_not_auto(params, "pcoupltype", args.pcoupltype)
    set_if_not_auto(params, "constraints", args.constraints)
    set_if_not_auto(params, "constraint-algorithm", args.constraint_algorithm)
    set_if_not_auto(params, "gen-vel", args.gen_vel)
    set_if_not_auto(params, "pull", args.pull)

    # Simple numeric/text keys (only if user provided a value)
    direct = [
        ("dt", args.dt),
        ("nsteps", args.nsteps),
        ("tinit", args.tinit),
        ("init-step", args.init_step),

        ("nstlist", args.nstlist),
        ("verlet-buffer-tolerance", args.verlet_buffer_tolerance),
        ("rlist", args.rlist),
        ("rcoulomb", args.rcoulomb),
        ("rvdw", args.rvdw),

        ("fourierspacing", args.fourierspacing),
        ("pme-order", args.pme_order),
        ("ewald-rtol", args.ewald_rtol),

        ("tc-grps", args.tc_grps),
        ("tau-t", args.tau_t),
        ("ref-t", args.ref_t),
        ("nsttcouple", args.nsttcouple),

        ("gen-temp", args.gen_temp),
        ("gen-seed", args.gen_seed),

        ("tau-p", args.tau_p),
        ("ref-p", args.ref_p),
        ("compressibility", args.compressibility),
        ("nstpcouple", args.nstpcouple),

        ("lincs-order", args.lincs_order),
        ("lincs-iter", args.lincs_iter),
        ("shake-tol", args.shake_tol),

        ("nstxout", args.nstxout),
        ("nstvout", args.nstvout),
        ("nstfout", args.nstfout),
        ("nstlog", args.nstlog),
        ("nstenergy", args.nstenergy),
        ("nstxout-compressed", args.nstxout_compressed),
        ("compressed-x-precision", args.compressed_x_precision),
    ]
    for k, v in direct:
        if v is not None and str(v).strip() != "":
            params[k] = str(v)

    # Pulling (basic coordinate 1)
    if args.pull_groups.strip():
        params["pull-ngroups"] = str(args.pull_ngroups).strip() if str(args.pull_ngroups).strip() else params.get("pull-ngroups", "1")
        params["pull-ncoords"] = str(args.pull_ncoords).strip() if str(args.pull_ncoords).strip() else params.get("pull-ncoords", "1")
        params["pull-coord1-groups"] = args.pull_groups.strip()

        if args.pull_coord1_type.strip():
            params["pull-coord1-type"] = args.pull_coord1_type.strip()
        if args.pull_coord1_geometry.strip():
            params["pull-coord1-geometry"] = args.pull_coord1_geometry.strip()
        if args.pull_coord1_dim.strip():
            params["pull-coord1-dim"] = args.pull_coord1_dim.strip()
        if args.pull_coord1_vec.strip():
            params["pull-coord1-vec"] = args.pull_coord1_vec.strip()
        if args.pull_coord1_init.strip():
            params["pull-coord1-init"] = args.pull_coord1_init.strip()
        if args.pull_coord1_rate.strip():
            params["pull-coord1-rate"] = args.pull_coord1_rate.strip()
        if args.pull_coord1_k.strip():
            params["pull-coord1-k"] = args.pull_coord1_k.strip()

    # Extra mdp lines override everything if key collisions
    extra_lines, extra_kv = parse_extra_lines(args.extra_mdp)
    for k, v in extra_kv.items():
        params[k] = v

    # Emit file
    with open(args.output, "w", encoding="utf-8") as f:
        f.write("; Generated by Galaxy tool: Generate GROMACS MDP File (Pro)\n")
        f.write(f"; Preset: {args.preset}\n")
        f.write(";\n")

        # Keep sections roughly organized for readability
        def write_section(title, keys):
            present = [k for k in keys if k in params]
            if not present:
                return
            f.write(f"\n; {title}\n")
            for k in present:
                f.write(f"{k:<22} = {params[k]}\n")

        write_section("Run control", ["integrator","tinit","dt","nsteps","init-step","continuation"])
        write_section("Neighbor searching", ["cutoff-scheme","nstlist","pbc","verlet-buffer-tolerance","rlist"])
        write_section("Electrostatics", ["coulombtype","coulomb-modifier","rcoulomb","fourierspacing","pme-order","ewald-rtol"])
        write_section("Van der Waals", ["vdwtype","vdw-modifier","rvdw","DispCorr"])
        write_section("Temperature coupling", ["tcoupl","tc-grps","tau-t","ref-t","nsttcouple"])
        write_section("Velocity generation", ["gen-vel","gen-temp","gen-seed"])
        write_section("Pressure coupling", ["pcoupl","pcoupltype","tau-p","ref-p","compressibility","nstpcouple"])
        write_section("Constraints", ["constraints","constraint-algorithm","lincs-order","lincs-iter","shake-tol"])
        write_section("Output control", ["nstxout","nstvout","nstfout","nstlog","nstenergy","nstxout-compressed","compressed-x-precision"])
        write_section("Pulling", ["pull","pull-ngroups","pull-ncoords","pull-coord1-groups","pull-coord1-type","pull-coord1-geometry",
                                 "pull-coord1-dim","pull-coord1-vec","pull-coord1-init","pull-coord1-rate","pull-coord1-k"])

        # Write remaining keys not in the above ordering, in insertion order
        ordered_keys = set([
            "integrator","tinit","dt","nsteps","init-step","continuation",
            "cutoff-scheme","nstlist","pbc","verlet-buffer-tolerance","rlist",
            "coulombtype","coulomb-modifier","rcoulomb","fourierspacing","pme-order","ewald-rtol",
            "vdwtype","vdw-modifier","rvdw","DispCorr",
            "tcoupl","tc-grps","tau-t","ref-t","nsttcouple",
            "gen-vel","gen-temp","gen-seed",
            "pcoupl","pcoupltype","tau-p","ref-p","compressibility","nstpcouple",
            "constraints","constraint-algorithm","lincs-order","lincs-iter","shake-tol",
            "nstxout","nstvout","nstfout","nstlog","nstenergy","nstxout-compressed","compressed-x-precision",
            "pull","pull-ngroups","pull-ncoords","pull-coord1-groups","pull-coord1-type","pull-coord1-geometry",
            "pull-coord1-dim","pull-coord1-vec","pull-coord1-init","pull-coord1-rate","pull-coord1-k"
        ])

        leftovers = [k for k in params.keys() if k not in ordered_keys]
        if leftovers:
            f.write("\n; Additional parameters\n")
            for k in leftovers:
                f.write(f"{k:<22} = {params[k]}\n")

        if args.extra_mdp.strip():
            f.write("\n; ---- Extra .mdp lines (verbatim as provided) ----\n")
            for line in extra_lines:
                f.write(line + "\n")

if __name__ == "__main__":
    main()

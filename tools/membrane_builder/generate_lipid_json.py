#!/usr/bin/env python3
"""
Script to parse lipid XML options and generate a comprehensive JSON file with lipid fragments.
"""

import re
import json

# All lipid definitions from macro.xml
lipid_data = """
        <option value="AHPA"> AHPA (1-arachidonoyl-2-docosahexaenoyl-sn-glycero-3-phosphate), Charge -1 </option>
        <option value="AHPC"> AHPC (1-arachidonoyl-2-docosahexaenoyl-sn-glycero-3-phosphocholine), Charge 0 </option>
        <option value="AHPE"> AHPE (1-arachidonoyl-2-docosahexaenoyl-sn-glycero-3-phosphoethanolamine), Charge 0 </option>
        <option value="AHPG"> AHPG (1-arachidonoyl-2-docosahexaenoyl-sn-glycero-3-phosphoglycerol), Charge -1 </option>
        <option value="AHPS"> AHPS (1-arachidonoyl-2-docosahexaenoyl-sn-glycero-3-phosphoserine), Charge -1 </option>
        <option value="ALPA"> ALPA (1-arachidonoyl-2-lauroyl-sn-glycero-3-phosphate), Charge -1 </option>
        <option value="ALPC"> ALPC (1-arachidonoyl-2-lauroyl-sn-glycero-3-phosphocholine), Charge 0 </option>
        <option value="ALPE"> ALPE (1-arachidonoyl-2-lauroyl-sn-glycero-3-phosphoethanolamine), Charge 0 </option>
        <option value="ALPG"> ALPG (1-arachidonoyl-2-lauroyl-sn-glycero-3-phosphoglycerol), Charge -1 </option>
        <option value="ALPS"> ALPS (1-arachidonoyl-2-lauroyl-sn-glycero-3-phosphoserine), Charge -1 </option>
        <option value="AMPA"> AMPA (1-arachidonoyl-2-myristoyl-sn-glycero-3-phosphate), Charge -1 </option>
        <option value="AMPC"> AMPC (1-arachidonoyl-2-myristoyl-sn-glycero-3-phosphocholine), Charge 0 </option>
        <option value="AMPE"> AMPE (1-arachidonoyl-2-myristoyl-sn-glycero-3-phosphoethanolamine), Charge 0 </option>
        <option value="AMPG"> AMPG (1-arachidonoyl-2-myristoyl-sn-glycero-3-phosphoglycerol), Charge -1 </option>
        <option value="AMPS"> AMPS (1-arachidonoyl-2-myristoyl-sn-glycero-3-phosphoserine), Charge -1 </option>
        <option value="AOPA"> AOPA (1-arachidonoyl-2-oleoyl-sn-glycero-3-phosphate), Charge -1 </option>
        <option value="AOPC"> AOPC (1-arachidonoyl-2-oleoyl-sn-glycero-3-phosphocholine), Charge 0 </option>
        <option value="AOPE"> AOPE (1-arachidonoyl-2-oleoyl-sn-glycero-3-phosphoethanolamine), Charge 0 </option>
        <option value="AOPG"> AOPG (1-arachidonoyl-2-oleoyl-sn-glycero-3-phosphoglycerol), Charge -1 </option>
        <option value="AOPS"> AOPS (1-arachidonoyl-2-oleoyl-sn-glycero-3-phosphoserine), Charge -1 </option>
        <option value="APPA"> APPA (1-arachidonoyl-2-palmitoyl-sn-glycero-3-phosphate), Charge -1 </option>
        <option value="APPC"> APPC (1-arachidonoyl-2-palmitoyl-sn-glycero-3-phosphocholine), Charge 0 </option>
        <option value="APPE"> APPE (1-arachidonoyl-2-palmitoyl-sn-glycero-3-phosphoethanolamine), Charge 0 </option>
        <option value="APPG"> APPG (1-arachidonoyl-2-palmitoyl-sn-glycero-3-phosphoglycerol), Charge -1 </option>
        <option value="APPS"> APPS (1-arachidonoyl-2-palmitoyl-sn-glycero-3-phosphoserine), Charge -1 </option>
        <option value="ASM"> ASM (N-arachidonoyl-D-erythro-sphingosylphosphorylcholine) [Only available with Lipid21], Charge 0 </option>
        <option value="ASPA"> ASPA (1-arachidonoyl-2-stearoyl-sn-glycero-3-phosphate), Charge -1 </option>
        <option value="ASPC"> ASPC (1-arachidonoyl-2-stearoyl-sn-glycero-3-phosphocholine), Charge 0 </option>
        <option value="ASPE"> ASPE (1-arachidonoyl-2-stearoyl-sn-glycero-3-phosphoethanolamine), Charge 0 </option>
        <option value="ASPG"> ASPG (1-arachidonoyl-2-stearoyl-sn-glycero-3-phosphoglycerol), Charge -1 </option>
        <option value="ASPS"> ASPS (1-arachidonoyl-2-stearoyl-sn-glycero-3-phosphoserine), Charge -1 </option>
        <option value="DAPA"> DAPA (1,2-diarachidonoyl-sn-glycero-3-phosphate), Charge -1 </option>
        <option value="DAPE"> DAPE (1,2-diarachidonoyl-sn-glycero-3-phosphoethanolamine), Charge 0 </option>
        <option value="DAPG"> DAPG (1,2-diarachidonoyl-sn-glycero-3-phosphoglycerol), Charge -1 </option>
        <option value="DAPS"> DAPS (1,2-diarachidonoyl-sn-glycero-3-phosphoserine), Charge -1 </option>
        <option value="DHPA"> DHPA (1,2-didocosahexaenoyl-sn-glycero-3-phosphate), Charge -1 </option>
        <option value="DHPC"> DHPC (1,2-didocosahexaenoyl-sn-glycero-3-phosphocholine), Charge 0 </option>
        <option value="DHPE"> DHPE (1,2-didocosahexaenoyl-sn-glycero-3-phosphoethanolamine), Charge 0 </option>
        <option value="DHPG"> DHPG (1,2-didocosahexaenoyl-sn-glycero-3-phosphoglycerol), Charge -1 </option>
        <option value="DHPS"> DHPS (1,2-didocosahexaenoyl-sn-glycero-3-phosphoserine), Charge -1 </option>
        <option value="DLPA"> DLPA (1,2-dilauroyl-sn-glycero-3-phosphate), Charge -1 </option>
        <option value="DLPC"> DLPC (1,2-dilauroyl-sn-glycero-3-phosphocholine), Charge 0 </option>
        <option value="DLPE"> DLPE (1,2-dilauroyl-sn-glycero-3-phosphoethanolamine), Charge 0 </option>
        <option value="DLPG"> DLPG (1,2-dilauroyl-sn-glycero-3-phosphoglycerol), Charge -1 </option>
        <option value="DLPS"> DLPS (1,2-dilauroyl-sn-glycero-3-phosphoserine), Charge -1 </option>
        <option value="DMPA"> DMPA (1,2-dimyristoyl-sn-glycero-3-phosphate), Charge -1 </option>
        <option value="DMPC"> DMPC (1,2-dimyristoyl-sn-glycero-3-phosphocholine), Charge 0 </option>
        <option value="DMPE"> DMPE (1,2-dimyristoyl-sn-glycero-3-phosphoethanolamine), Charge 0 </option>
        <option value="DMPG"> DMPG (1,2-dimyristoyl-sn-glycero-3-phosphoglycerol), Charge -1 </option>
        <option value="DMPS"> DMPS (1,2-dimyristoyl-sn-glycero-3-phosphoserine), Charge -1 </option>
        <option value="DOPA"> DOPA (1,2-dioleoyl-sn-glycero-3-phosphate), Charge -1 </option>
        <option value="DOPC"> DOPC (1,2-dioleoyl-sn-glycero-3-phosphocholine), Charge 0 </option>
        <option value="DOPE"> DOPE (1,2-dioleoyl-sn-glycero-3-phosphoethanolamine), Charge 0 </option>
        <option value="DOPG"> DOPG (1,2-dioleoyl-sn-glycero-3-phosphoglycerol), Charge -1 </option>
        <option value="DOPS"> DOPS (1,2-dioleoyl-sn-glycero-3-phosphoserine), Charge -1 </option>
        <option value="DPPA"> DPPA (1,2-dipalmitoyl-sn-glycero-3-phosphate), Charge -1 </option>
        <option value="DPPE"> DPPE (1,2-dipalmitoyl-sn-glycero-3-phosphoethanolamine), Charge 0 </option>
        <option value="DPPG"> DPPG (1,2-dipalmitoyl-sn-glycero-3-phosphoglycerol), Charge -1 </option>
        <option value="DPPS"> DPPS (1,2-dipalmitoyl-sn-glycero-3-phosphoserine), Charge -1 </option>
        <option value="DSPA"> DSPA (1,2-distearoyl-sn-glycero-3-phosphate), Charge -1 </option>
        <option value="DSPC"> DSPC (1,2-distearoyl-sn-glycero-3-phosphocholine), Charge 0 </option>
        <option value="DSPE"> DSPE (1,2-distearoyl-sn-glycero-3-phosphoethanolamine), Charge 0 </option>
        <option value="DSPG"> DSPG (1,2-distearoyl-sn-glycero-3-phosphoglycerol), Charge -1 </option>
        <option value="DSPS"> DSPS (1,2-distearoyl-sn-glycero-3-phosphoserine), Charge -1 </option>
        <option value="HAPA"> HAPA (1-docosahexaenoyl-2-arachidonoyl-sn-glycero-3-phosphate), Charge -1 </option>
        <option value="HAPC"> HAPC (1-docosahexaenoyl-2-arachidonoyl-sn-glycero-3-phosphocholine), Charge 0 </option>
        <option value="HAPE"> HAPE (1-docosahexaenoyl-2-arachidonoyl-sn-glycero-3-phosphoethanolamine), Charge 0 </option>
        <option value="HAPG"> HAPG (1-docosahexaenoyl-2-arachidonoyl-sn-glycero-3-phosphoglycerol), Charge -1 </option>
        <option value="HAPS"> HAPS (1-docosahexaenoyl-2-arachidonoyl-sn-glycero-3-phosphoserine), Charge -1 </option>
        <option value="HLPA"> HLPA (1-docosahexaenoyl-2-lauroyl-sn-glycero-3-phosphate), Charge -1 </option>
        <option value="HLPC"> HLPC (1-docosahexaenoyl-2-lauroyl-sn-glycero-3-phosphocholine), Charge 0 </option>
        <option value="HLPE"> HLPE (1-docosahexaenoyl-2-lauroyl-sn-glycero-3-phosphoethanolamine), Charge 0 </option>
        <option value="HLPG"> HLPG (1-docosahexaenoyl-2-lauroyl-sn-glycero-3-phosphoglycerol), Charge -1 </option>
        <option value="HLPS"> HLPS (1-docosahexaenoyl-2-lauroyl-sn-glycero-3-phosphoserine), Charge -1 </option>
        <option value="HMPA"> HMPA (1-docosahexaenoyl-2-myristoyl-sn-glycero-3-phosphate), Charge -1 </option>
        <option value="HMPC"> HMPC (1-docosahexaenoyl-2-myristoyl-sn-glycero-3-phosphocholine), Charge 0 </option>
        <option value="HMPE"> HMPE (1-docosahexaenoyl-2-myristoyl-sn-glycero-3-phosphoethanolamine), Charge 0 </option>
        <option value="HMPG"> HMPG (1-docosahexaenoyl-2-myristoyl-sn-glycero-3-phosphoglycerol), Charge -1 </option>
        <option value="HMPS"> HMPS (1-docosahexaenoyl-2-myristoyl-sn-glycero-3-phosphoserine), Charge -1 </option>
        <option value="HOPA"> HOPA (1-docosahexaenoyl-2-oleoyl-sn-glycero-3-phosphate), Charge -1 </option>
        <option value="HOPC"> HOPC (1-docosahexaenoyl-2-oleoyl-sn-glycero-3-phosphocholine), Charge 0 </option>
        <option value="HOPE"> HOPE (1-docosahexaenoyl-2-oleoyl-sn-glycero-3-phosphoethanolamine), Charge 0 </option>
        <option value="HOPG"> HOPG (1-docosahexaenoyl-2-oleoyl-sn-glycero-3-phosphoglycerol), Charge -1 </option>
        <option value="HOPS"> HOPS (1-docosahexaenoyl-2-oleoyl-sn-glycero-3-phosphoserine), Charge -1 </option>
        <option value="HPPA"> HPPA (1-docosahexaenoyl-2-palmitoyl-sn-glycero-3-phosphate), Charge -1 </option>
        <option value="HPPC"> HPPC (1-docosahexaenoyl-2-palmitoyl-sn-glycero-3-phosphocholine), Charge 0 </option>
        <option value="HPPE"> HPPE (1-docosahexaenoyl-2-palmitoyl-sn-glycero-3-phosphoethanolamine), Charge 0 </option>
        <option value="HPPG"> HPPG (1-docosahexaenoyl-2-palmitoyl-sn-glycero-3-phosphoglycerol), Charge -1 </option>
        <option value="HPPS"> HPPS (1-docosahexaenoyl-2-palmitoyl-sn-glycero-3-phosphoserine), Charge -1 </option>
        <option value="HSM"> HSM (N-docosahexaenoyl-D-erythro-sphingosylphosphorylcholine) [Only available with Lipid21], Charge 0 </option>
        <option value="HSPA"> HSPA (1-docosahexaenoyl-2-stearoyl-sn-glycero-3-phosphate), Charge -1 </option>
        <option value="HSPC"> HSPC (1-docosahexaenoyl-2-stearoyl-sn-glycero-3-phosphocholine), Charge 0 </option>
        <option value="HSPE"> HSPE (1-docosahexaenoyl-2-stearoyl-sn-glycero-3-phosphoethanolamine), Charge 0 </option>
        <option value="HSPG"> HSPG (1-docosahexaenoyl-2-stearoyl-sn-glycero-3-phosphoglycerol), Charge -1 </option>
        <option value="HSPS"> HSPS (1-docosahexaenoyl-2-stearoyl-sn-glycero-3-phosphoserine), Charge -1 </option>
        <option value="LAPA"> LAPA (1-lauroyl-2-arachidonoyl-sn-glycero-3-phosphate), Charge -1 </option>
        <option value="LAPC"> LAPC (1-lauroyl-2-arachidonoyl-sn-glycero-3-phosphocholine), Charge 0 </option>
        <option value="LAPE"> LAPE (1-lauroyl-2-arachidonoyl-sn-glycero-3-phosphoethanolamine), Charge 0 </option>
        <option value="LAPG"> LAPG (1-lauroyl-2-arachidonoyl-sn-glycero-3-phosphoglycerol), Charge -1 </option>
        <option value="LAPS"> LAPS (1-lauroyl-2-arachidonoyl-sn-glycero-3-phosphoserine), Charge -1 </option>
        <option value="LHPA"> LHPA (1-lauroyl-2-docosahexaenoyl-sn-glycero-3-phosphate), Charge -1 </option>
        <option value="LHPC"> LHPC (1-lauroyl-2-docosahexaenoyl-sn-glycero-3-phosphocholine), Charge 0 </option>
        <option value="LHPE"> LHPE (1-lauroyl-2-docosahexaenoyl-sn-glycero-3-phosphoethanolamine), Charge 0 </option>
        <option value="LHPG"> LHPG (1-lauroyl-2-docosahexaenoyl-sn-glycero-3-phosphoglycerol), Charge -1 </option>
        <option value="LHPS"> LHPS (1-lauroyl-2-docosahexaenoyl-sn-glycero-3-phosphoserine), Charge -1 </option>
        <option value="LMPA"> LMPA (1-lauroyl-2-myristoyl-sn-glycero-3-phosphate), Charge -1 </option>
        <option value="LMPC"> LMPC (1-lauroyl-2-myristoyl-sn-glycero-3-phosphocholine), Charge 0 </option>
        <option value="LMPE"> LMPE (1-lauroyl-2-myristoyl-sn-glycero-3-phosphoethanolamine), Charge 0 </option>
        <option value="LMPG"> LMPG (1-lauroyl-2-myristoyl-sn-glycero-3-phosphoglycerol), Charge -1 </option>
        <option value="LMPS"> LMPS (1-lauroyl-2-myristoyl-sn-glycero-3-phosphoserine), Charge -1 </option>
        <option value="LOPA"> LOPA (1-lauroyl-2-oleoyl-sn-glycero-3-phosphate), Charge -1 </option>
        <option value="LOPC"> LOPC (1-lauroyl-2-oleoyl-sn-glycero-3-phosphocholine), Charge 0 </option>
        <option value="LOPE"> LOPE (1-lauroyl-2-oleoyl-sn-glycero-3-phosphoethanolamine), Charge 0 </option>
        <option value="LOPG"> LOPG (1-lauroyl-2-oleoyl-sn-glycero-3-phosphoglycerol), Charge -1 </option>
        <option value="LOPS"> LOPS (1-lauroyl-2-oleoyl-sn-glycero-3-phosphoserine), Charge -1 </option>
        <option value="LPPA"> LPPA (1-lauroyl-2-palmitoyl-sn-glycero-3-phosphate), Charge -1 </option>
        <option value="LPPC"> LPPC (1-lauroyl-2-palmitoyl-sn-glycero-3-phosphocholine), Charge 0 </option>
        <option value="LPPE"> LPPE (1-lauroyl-2-palmitoyl-sn-glycero-3-phosphoethanolamine), Charge 0 </option>
        <option value="LPPG"> LPPG (1-lauroyl-2-palmitoyl-sn-glycero-3-phosphoglycerol), Charge -1 </option>
        <option value="LPPS"> LPPS (1-lauroyl-2-palmitoyl-sn-glycero-3-phosphoserine), Charge -1 </option>
        <option value="LSM"> LSM (N-lauroyl-D-erythro-sphingosylphosphorylcholine) [Only available with Lipid21], Charge 0 </option>
        <option value="LSPA"> LSPA (1-lauroyl-2-stearoyl-sn-glycero-3-phosphate), Charge -1 </option>
        <option value="LSPC"> LSPC (1-lauroyl-2-stearoyl-sn-glycero-3-phosphocholine), Charge 0 </option>
        <option value="LSPE"> LSPE (1-lauroyl-2-stearoyl-sn-glycero-3-phosphoethanolamine), Charge 0 </option>
        <option value="LSPG"> LSPG (1-lauroyl-2-stearoyl-sn-glycero-3-phosphoglycerol), Charge -1 </option>
        <option value="LSPS"> LSPS (1-lauroyl-2-stearoyl-sn-glycero-3-phosphoserine), Charge -1 </option>
        <option value="POPC"> POPC (1-palmitoyl-2-oleoyl-sn-glycero-3-phosphocholine), Charge 0 </option>
        <option value="POPE"> POPE (1-palmitoyl-2-oleoyl-sn-glycero-3-phosphoethanolamine), Charge 0 </option>
        <option value="POPG"> POPG (1-palmitoyl-2-oleoyl-sn-glycero-3-phosphoglycerol), Charge -1 </option>
        <option value="POPS"> POPS (1-palmitoyl-2-oleoyl-sn-glycero-3-phosphoserine), Charge -1 </option>
        <option value="PSM"> PSM (N-palmitoyl-D-erythro-sphingosylphosphorylcholine) [Only available with Lipid21], Charge 0 </option>
        <option value="siDMPC"> siDMPC (1,2-dimyristoyl-sn-glycero-3-phosphocholine) [Only available with SIRAH], Charge 0 </option>
        <option value="siDMPE"> siDMPE (1,2-dimyristoyl-sn-glycero-3-phosphoethanolamine) [Only available with SIRAH], Charge 0 </option>
        <option value="siDMPS"> siDMPS (1,2-dimyristoyl-sn-glycero-3-phosphoserine) [Only available with SIRAH], Charge -1 </option>
        <option value="siDOPC"> siDOPC (1,2-dioleoyl-sn-glycero-3-phosphocholine) [Only available with SIRAH], Charge 0 </option>
        <option value="siDOPE"> siDOPE (1,2-dioleoyl-sn-glycero-3-phosphoethanolamine) [Only available with SIRAH], Charge 0 </option>
        <option value="siDOPS"> siDOPS (1,2-dioleoyl-sn-glycero-3-phosphoserine) [Only available with SIRAH], Charge -1 </option>
        <option value="siDPPC"> siDPPC (1,2-dipalmitoyl-sn-glycero-3-phosphocholine) [Only available with SIRAH], Charge 0 </option>
        <option value="siDPPE"> siDPPE (1,2-dipalmitoyl-sn-glycero-3-phosphoethanolamine) [Only available with SIRAH], Charge 0 </option>
        <option value="siDPPS"> siDPPS (1,2-dipalmitoyl-sn-glycero-3-phosphoserine) [Only available with SIRAH], Charge -1 </option>
"""

# Acyl chain abbreviations
acyl_chains = {
    'A': 'arachidonoyl',
    'H': 'docosahexaenoyl',
    'L': 'lauroyl',
    'M': 'myristoyl',
    'O': 'oleoyl',
    'P': 'palmitoyl',
    'S': 'stearoyl',
}

# Head group abbreviations
head_groups = {
    'PA': ('phosphate', 'PA'),
    'PC': ('phosphocholine', 'PC'),
    'PE': ('phosphoethanolamine', 'PE'),
    'PG': ('phosphoglycerol', 'PG'),
    'PS': ('phosphoserine', 'PS'),
    'SM': ('sphingomyelin', 'SM'),
}

def extract_fragments(lipid_code):
    """Extract fragments from lipid code."""
    # Handle special cases like ASM, LSM, PSM, HSM
    if lipid_code.endswith('SM'):
        acyl = lipid_code[:-2]
        return [acyl, 'SM']
    
    # Handle Di- lipids (DX where X is acyl type)
    if lipid_code.startswith('D'):
        acyl = 'D' + lipid_code[1]
        head = lipid_code[2:]
        return [acyl, head]
    
    # Handle si- lipids (SIRAH)
    if lipid_code.startswith('si'):
        # siDMPC -> si + DM + PC
        rest = lipid_code[2:]
        if rest.startswith('D'):
            acyl = 'D' + rest[1]
            head = rest[2:]
            return ['si' + acyl, head]
    
    # Standard two acyl chain lipids
    if len(lipid_code) >= 3:
        acyl1 = lipid_code[0]
        acyl2 = lipid_code[1]
        head = lipid_code[2:]
        return [acyl1 + acyl2, head]
    
    return []

def parse_lipid_line(line):
    """Parse a single lipid option line."""
    match = re.search(r'value="([^"]+)">\s*\1\s*\(([^)]+)\).*?Charge\s*([-\d]+)', line)
    if match:
        lipid_code = match.group(1)
        description = match.group(2)
        charge = int(match.group(3))
        
        fragments = extract_fragments(lipid_code)
        
        return {
            'name': lipid_code,
            'description': description,
            'charge': charge,
            'fragments': fragments
        }
    return None

def main():
    """Parse all lipids and create JSON."""
    lipids = {}
    
    for line in lipid_data.split('\n'):
        if '<option value=' in line:
            lipid_info = parse_lipid_line(line)
            if lipid_info:
                lipids[lipid_info['name']] = {
                    'name': lipid_info['name'],
                    'description': lipid_info['description'],
                    'charge': lipid_info['charge'],
                    'fragments': lipid_info['fragments']
                }
    
    # Write to JSON file
    output_file = '/home/joshij/Desktop/MemGen/membrane_builder/lipids_complete.json'
    with open(output_file, 'w') as f:
        json.dump(lipids, f, indent=2)
    
    print(f"Successfully created {output_file}")
    print(f"Total lipids parsed: {len(lipids)}")
    
    # Print some examples
    print("\nFirst 5 lipids:")
    for i, (key, value) in enumerate(list(lipids.items())[:5]):
        print(f"  {key}: fragments={value['fragments']}, charge={value['charge']}")

if __name__ == '__main__':
    main()

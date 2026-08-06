"""
Sort planpapers/ into folders by project direction, and give the files readable
names.

Folders mirror LITERATURE.md. Where a paper could sit in two places it goes where
it is most likely to be LOOKED FOR, not where it is most cited:

  * Takano 2014 uses ambient noise but reports a tidal result -> 02_tides
  * Rubinstein & Beroza 2005 uses repeaters but is a Parkfield site result -> 05_site
  * Lellouch's two 2019 papers are DAS method papers, but their value here is that
    they describe THIS cable -> 05_site
  * Ben-Zion & Leary and Berger are thermoelastic, not tidal, but thermoelastic
    strain is the confound a tidal measurement must exclude -> 02_tides

Renaming to author_year_journal_topic.pdf is the other half. Names like
"science.285.5428.718.pdf" hide the fact that the file is Nadeau & McEvilly 1999.
Every move is recorded in MANIFEST.csv, so the whole operation reverses with a
two-line script.

Run with DRY=1 to preview.
"""
import csv
import os
import shutil
import unicodedata

ROOT = '/home/groups/ettore88/nberrios/planpapers'
DRY = os.environ.get('DRY', '0') == '1'

FOLDERS = {
    '01_repeaters': 'Repeating earthquakes -> fault slip at depth',
    '02_tides': 'Tides -> velocity modulation and poroelastic response',
    '03_ambient_noise': 'Ambient-noise cross-correlation -> continuous monitoring',
    '04_methods_das': 'DAS as an instrument (cross-cutting method)',
    '05_site_safod_parkfield': 'SAFOD / Parkfield -- what is known about this hole',
    '06_methods_borehole': 'Borehole and active-source methods',
    '_duplicates': 'Second copies, kept rather than deleted',
    '_not_papers': 'Not literature',
}

# old filename (without .pdf unless noted) -> (folder, new name)
MAP = {
    # ---------------- 01 repeaters ----------------
    'science.285.5428.718':
        ('01_repeaters', 'nadeau_mcevilly_1999_science_fault-slip-rates-at-depth'),
    'bssa0880030790':
        ('01_repeaters', 'nadeau_johnson_1998_bssa_parkfield-VI-moment-release-scaling'),
    'science.267.5197.503':
        ('01_repeaters', 'nadeau_1995_science_clustering-periodic-recurrence-parkfield'),
    'Journal of Geophysical Research  Solid Earth - 10 July 1984 - Poupinet - Monitoring velocity variations in the crust using':
        ('01_repeaters', 'poupinet_ellsworth_frechet_1984_jgr_earthquake-doublets-calaveras'),
    'annurev-earth-053018-060119':
        ('01_repeaters', 'uchida_burgmann_2019_annurev_repeating-earthquakes-review'),
    's40645-019-0284-z':
        ('01_repeaters', 'uchida_2019_peps_detection-of-repeating-earthquakes'),
    'Journal of Geophysical Research  Solid Earth - 2001 - Beeler - Earthquake stress drop and laboratory‐inferred interseismic':
        ('01_repeaters', 'beeler_hickman_2001_jgr_stress-drop-strength-recovery'),
    'Journal of Geophysical Research  Solid Earth - 2009 - Chen - Scaling of small repeating earthquakes explained by':
        ('01_repeaters', 'chen_lapusta_2009_jgr_scaling-repeaters-rate-and-state'),
    'Geophysical Research Letters - 2014 - Abercrombie - Stress drops of repeating earthquakes on the San Andreas Fault at':
        ('01_repeaters', 'abercrombie_2014_grl_stress-drops-repeating-earthquakes'),
    'Geophysical Research Letters - 2021 - Gao - Misconception of Waveform Similarity in the Identification of Repeating':
        ('01_repeaters', 'gao_2021_grl_misconception-waveform-similarity'),
    'Journal of Geophysical Research  Solid Earth - 2009 - Lengliné - Inferring the coseismic and postseismic stress changes':
        ('01_repeaters', 'lengline_marsan_2009_jgr_stress-changes-parkfield-repeaters'),
    'nature05780':
        ('01_repeaters', 'ide_2007_nature_scaling-law-for-slow-earthquakes'),

    # ---------------- 02 tides ----------------
    'Journal of Geophysical Research  1896-1977 - 10 March 1973 - De Fazio - Solid earth tide and observed change in the in situ':
        ('02_tides', 'defazio_aki_alba_1973_jgr_solid-earth-tide-velocity-change'),
    'Geophysical Research Letters - 2014 - Takano - Seismic velocity changes caused by the Earth tide  Ambient noise correlation (1)':
        ('02_tides', 'takano_2014_grl_velocity-changes-earth-tide-ambient-noise'),
    'Journal of Geophysical Research  1896-1977 - 10 January 1975 - Berger - A note on thermoelastic strains and tilts':
        ('02_tides', 'berger_1975_jgr_thermoelastic-strains-and-tilts'),
    'bssa0760051447':
        ('02_tides', 'benzion_leary_1986_bssa_thermoelastic-strain-halfspace'),
    'Journal of Geophysical Research  Solid Earth - 2008 - Métivier - Body tides of a convecting  laterally heterogeneous  and':
        ('02_tides', 'metivier_conrad_2008_jgr_body-tides-heterogeneous-earth'),
    'Journal of Geophysical Research  Solid Earth - 10 September 1989 - Rojstaczer - The influence of formation material':
        ('02_tides', 'rojstaczer_1989_jgr_well-water-levels-response-earth-tides'),
    'hess-26-4301-2022':
        ('02_tides', 'hess_2022_groundwater-response-earth-atmospheric-tides'),

    # ---------------- 03 ambient noise ----------------
    'AmbientNoiseTutorial':
        ('03_ambient_noise', 'wapenaar_2010_geophysics_interferometry-tutorial-part1'),
    'Geophysical Research Letters - 2004 - Shapiro - Emergence of broadband Rayleigh waves from correlations of the ambient':
        ('03_ambient_noise', 'shapiro_campillo_2004_grl_broadband-rayleigh-from-noise'),
    'science.1160943':
        ('03_ambient_noise', 'brenguier_2008_science_postseismic-relaxation-parkfield'),
    'JGR Solid Earth - 2023 - Li - Daily and Seasonal Variations of Shallow Seismic Velocities in Southern California From Joint':
        ('03_ambient_noise', 'li_benzion_2023_jgr_daily-seasonal-shallow-velocity'),

    # ---------------- 04 DAS as instrument ----------------
    'srl-2019112.1':
        ('04_methods_das', 'zhan_2020_srl_das-fiber-optic-seismic-antennas'),
    'JGR Solid Earth - 2020 - Lindsey - On the Broadband Instrument Response of Fiber‐Optic DAS Arrays':
        ('04_methods_das', 'lindsey_2020_jgr_broadband-instrument-response-das'),
    'JGR Solid Earth - 2022 - Ichinose - Comparisons Between Array Derived Dynamic Strain Rate  ADDS  and Fiber‐Optic':
        ('04_methods_das', 'ichinose_2022_jgr_das-vs-array-derived-strain-rate'),
    'srl-2020149.1':
        ('04_methods_das', 'lellouch_2020_srl_das-vs-geophones-forge'),
    'tle35070610.1':
        ('04_methods_das', 'madsen_2016_tle_data-driven-depth-calibration-das'),
    'dissertation':
        ('04_methods_das', 'martin_2018_thesis_passive-imaging-das'),
    'JGR Solid Earth - 2024 - Atterholt - Imaging the Garlock Fault Zone With a Fiber  A Limited Damage Zone and Hidden':
        ('04_methods_das', 'atterholt_2024_jgr_garlock-fault-zone-with-fiber'),

    # ---------------- 05 SAFOD / Parkfield site ----------------
    'JGR Solid Earth - 2019 - Lellouch - Seismic Velocity Estimation Using Passive Downhole Distributed Acoustic Sensing Records':
        ('05_site_safod_parkfield', 'lellouch_2019_jgr_velocity-estimation-downhole-das-SAFOD'),
    'bssa-2019176.1 (2)':
        ('05_site_safod_parkfield', 'lellouch_2019_bssa_velocity-based-detection-SAFOD'),
    'science.1090711':
        ('05_site_safod_parkfield', 'chavarria_2003_science_vsp-inside-san-andreas-parkfield'),
    'Geophysical Research Letters - 2004 - Hickman - Stress orientations and magnitudes in the SAFOD pilot hole':
        ('05_site_safod_parkfield', 'hickman_zoback_2004_grl_stress-orientations-pilot-hole'),
    'Geophysical Research Letters - 2004 - Boness - Stress‐induced seismic velocity anisotropy and physical properties in the':
        ('05_site_safod_parkfield', 'boness_zoback_2004_grl_stress-induced-velocity-anisotropy'),
    'Geophysical Research Letters - 2004 - Townend - Regional tectonic stress near the San Andreas fault in central and southern':
        ('05_site_safod_parkfield', 'townend_zoback_2004_grl_regional-tectonic-stress'),
    'Geophysical Research Letters - 2004 - Li - Low‐velocity damaged structure of the San Andreas Fault at Parkfield from fault':
        ('05_site_safod_parkfield', 'li_2004_grl_low-velocity-damaged-structure-trapped-waves'),
    'Geophysical Research Letters - 2004 - Unsworth - Electrical resistivity structure at the SAFOD site from magnetotelluric':
        ('05_site_safod_parkfield', 'unsworth_2004_grl_magnetotelluric-resistivity-SAFOD'),
    'Geophysical Research Letters - 2006 - Hole - Structure of the San Andreas fault zone at SAFOD from a seismic refraction':
        ('05_site_safod_parkfield', 'hole_2006_grl_fault-zone-structure-seismic-refraction'),
    'Geophysical Research Letters - 2006 - Schleicher - Origin and significance of clay‐coated fractures in mudrock fragments of':
        ('05_site_safod_parkfield', 'schleicher_2006_grl_clay-coated-fractures-SAFOD'),
    'Geophysical Research Letters - 2005 - Rubinstein - Depth constraints on nonlinear strong ground motion from the 2004':
        ('05_site_safod_parkfield', 'rubinstein_beroza_2005_grl_depth-constraints-nonlinear-ground-motion'),
    'Ellsworth_and Malin_Sibson_Volume':
        ('05_site_safod_parkfield', 'ellsworth_malin_geolsoc_deep-rock-damage-guided-waves'),
    'Cochran_2025':
        ('05_site_safod_parkfield', 'cochran_2025_seismica_continental-scientific-drilling'),
    'DrillbitseismicimagesfracturesofSanAndreasfaultsystem-OilGasJournal':
        ('05_site_safod_parkfield', 'oilgasjournal_2005_drillbit-seismic-san-andreas_TRADE-ARTICLE'),

    # ---------------- 06 borehole / active source ----------------
    'Near Surface Geophysics - 2021 - Banerjee - Fracture analysis using Stoneley waves in a coalbed methane reservoir':
        ('06_methods_borehole', 'banerjee_chatterjee_2021_nsg_stoneley-waves-fracture-analysis'),

    # ---------------- housekeeping ----------------
    'science.1090711 (1)':
        ('_duplicates', 'chavarria_2003_science_vsp-inside-san-andreas_COPY2'),
}

NOT_PAPERS = {'awd_borehole_coupling.ipynb': ('_not_papers',
                                              'awd_borehole_coupling.ipynb')}


def _n(s):
    """NFC-normalise. Accented filenames (Metivier, Lengline) are stored
    decomposed on this filesystem but composed in the map above, so a plain
    string compare misses them."""
    return unicodedata.normalize('NFC', s)


def main():
    # disk name keyed by normalised form, so lookups survive the NFC/NFD split
    disk = {_n(f[:-4]): f[:-4] for f in os.listdir(ROOT) if f.endswith('.pdf')}
    present = set(disk)
    mapped = {_n(k) for k in MAP}
    missing = mapped - present
    unmapped = present - mapped

    if missing:
        print(f'IN MAP BUT NOT ON DISK ({len(missing)}):')
        for m in sorted(missing):
            print(f'   {m}')
    if unmapped:
        print(f'\nON DISK BUT NOT MAPPED ({len(unmapped)}) -- these would be left '
              f'in place:')
        for m in sorted(unmapped):
            print(f'   {m}')
    if missing or unmapped:
        print('\nresolve the above before moving; nothing has been touched.')
        if missing:
            return

    for d, desc in FOLDERS.items():
        p = os.path.join(ROOT, d)
        if not DRY:
            os.makedirs(p, exist_ok=True)
            with open(os.path.join(p, 'README.md'), 'w') as fh:
                fh.write(f'# {d}\n\n{desc}\n\n'
                         f'See `notebooks/LITERATURE.md` in the safod-das repo for '
                         f'why each paper is here and what it is for.\n')

    rows = []
    for old, (folder, new) in sorted(MAP.items(), key=lambda x: (x[1][0], x[1][1])):
        old = disk.get(_n(old), old)          # use the on-disk spelling
        src = os.path.join(ROOT, old + '.pdf')
        dst = os.path.join(ROOT, folder, new + '.pdf')
        if not os.path.exists(src):
            continue
        print(f'{folder:26s} {new}.pdf')
        if not DRY:
            shutil.move(src, dst)
        rows.append(dict(folder=folder, new_name=new + '.pdf',
                         old_name=old + '.pdf'))

    for old, (folder, new) in NOT_PAPERS.items():
        src = os.path.join(ROOT, old)
        if os.path.exists(src):
            print(f'{folder:26s} {new}')
            if not DRY:
                shutil.move(src, os.path.join(ROOT, folder, new))
            rows.append(dict(folder=folder, new_name=new, old_name=old))

    if not DRY:
        with open(os.path.join(ROOT, 'MANIFEST.csv'), 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=['folder', 'new_name', 'old_name'])
            w.writeheader()
            w.writerows(rows)

    print(f'\n{len(rows)} files {"would move" if DRY else "moved"}')
    for d in FOLDERS:
        p = os.path.join(ROOT, d)
        n = len([x for x in os.listdir(p) if x.endswith(('.pdf', '.ipynb'))]) \
            if os.path.isdir(p) else 0
        print(f'  {d:26s} {n:2d}')
    if not DRY:
        print(f'\nMANIFEST.csv records every move. To reverse:')
        print(f'  python3 -c "import csv,shutil,os;R=\'{ROOT}\';'
              f'[shutil.move(os.path.join(R,r[\'folder\'],r[\'new_name\']),'
              f'os.path.join(R,r[\'old_name\'])) '
              f'for r in csv.DictReader(open(os.path.join(R,\'MANIFEST.csv\')))]"')


if __name__ == '__main__':
    main()

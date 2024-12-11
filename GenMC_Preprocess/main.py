import yaml
import compile_vasp_data as cvd
import mag_distrib as mag


with open('param_in', 'r') as filehandle:
    param = yaml.safe_load(filehandle)

do_compile_data = param['do_compile_data']
root_dir = param['root_dir']
data_out = param['write_data']
species = param['species']
read_mag = param['read_mag']
use_spin_tol = param['use_spin_tol']
spin_tol = param['spin_tol']

do_mag_distrib = param['do_mag_distrib']
mag_out = param['write_mag_distrib']

do_mag_pair = param['do_mag_pair']
mag_pair_out = param['write_mag_pair']
target = param['target']
pair = param['pair']
clust = param['clust']

if do_compile_data:
    cvd.compile_vasp(root_dir, data_out, species, read_mag, use_spin_tol, spin_tol)

if do_mag_distrib:
    mag.find_mag_distrib(root_dir, mag_out, species)

if do_mag_pair:
    mag.find_pair_mag(root_dir, mag_pair_out, target, pair, clust)
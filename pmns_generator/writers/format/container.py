import pickle
from pathlib import Path
import sys
from sage.all import ceil, block_matrix, matrix

# Gestion des imports relatifs/système
CURRENT_DIR = Path(__file__).resolve().parent
ROOT_DIR = CURRENT_DIR.parent.parent.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.append(str(ROOT_DIR))
    
OUTPUT_DIR = "pmns_exec"
OUTPUT_DIR_SAVES = Path(OUTPUT_DIR) / "saves"

from ..format import format_element as format

from pmns_factory.core.operations.convertions_gestion import compute_conversion_tables, compute_nb_internal_reductions, gen_transition_matrix
from pmns_factory.core.operations.reductions.babai_reduction import gen_params_for_babai
from pmns_factory.core.operations.reductions.montgomery_reduction import gen_mn_reduction_matrix, search_polynomial_m
from pmns_factory.core.parameters.matrix_gestion import gen_overflow_matrix

class PMNSContainer:
    def __init__(self, Etype:int, pmns: dict):
        self.params = {**pmns}
        self.params['type'] = Etype
        self._build()



    def _ensure_list(self, obj):
        if hasattr(obj, 'rows'):
            return [list(row) for row in obj.rows()]
        return obj



    def _build(self):
        E, L = self.params['E'], self.params['L']
        p, k = self.params['p'], self.params['k']
        phi_pow = self.params['phi_pow']
        phi = 2**phi_pow
        gamma = self.params['gamma']
        rho = self.params['rho']
        n = E.degree()
        
        self.params['n'] = n
        self.params['n_limbs'] = ceil(p.nbits() / phi_pow)
        
        self.params['L'] = self._ensure_list(L)
        self.params['L_origin'] = L
        L_inv = -L.inverse() % phi
        self.params['L_inv'] = self._ensure_list(L_inv)
        self.params['L_inv_origin'] = L_inv

        ext_mat_sage = block_matrix([[gen_overflow_matrix(E)], [matrix(1, n)]])
        self.params['ext_red_mat'] = self._ensure_list(ext_mat_sage)
        
        M = search_polynomial_m(L, k, p, gamma, E)
        M_mat, N_mat = gen_mn_reduction_matrix(M, E, phi)
        self.params['M_mat'] = self._ensure_list(M_mat)
        self.params['N_mat'] = self._ensure_list(N_mat)
        
        h1, h2, L_inv_babai = gen_params_for_babai(L, phi_pow, rho, E)
        self.params['h1'] = h1
        self.params['h2'] = h2
        self.params['L_inv_babai'] = self._ensure_list(L_inv_babai)
        self.params['L_inv_babai_origin'] = L_inv_babai
        
        T_mat = gen_transition_matrix(gamma, k)
        self.params['T_mat'] = self._ensure_list(T_mat)
        self.params['T_mat_origin'] = T_mat
        
        
        
        
    def add_conversions_parameters(self):
        n, L = self.params['n'], self.params['L_origin']
        p, k = self.params['p'], self.params['k']
        phi_pow = self.params['phi_pow']
        phi = 2**phi_pow
        rho = self.params['rho']
        
        theta_pow = ceil(p.nbits() * k / n)
        n_red_extact = compute_nb_internal_reductions((2*rho)**(n/k), phi, rho, L)
        n_red_pseudo = compute_nb_internal_reductions(n * 2**ceil(p.nbits() * k / n) * (1/2 *L.norm(1))**2 , phi, rho, L)
        
        n_pol = ceil(n/k)
        self.params["theta_pow"] = theta_pow     
        self.params['n_pol'] = n_pol
        
        self.params['n_red_exact'] = n_red_extact
        self.params['n_red_pseudo'] = n_red_pseudo
        self.params['n_red_fast'] = 1
        
        self.params['fast_pols'] = compute_conversion_tables(self, theta_pow, 1, n_pol, over_field=True)
        self.params['int_pols'] = compute_conversion_tables(self, theta_pow, n_red_pseudo, n_pol, over_field=False)
        z_pols = compute_conversion_tables(self, 0, n_red_pseudo, 1, over_field=True)
        self.params['z_pols'] = [row[0] for row in z_pols]



    def save(self):
        path = OUTPUT_DIR_SAVES
        path.mkdir(parents=True, exist_ok=True)
        
        filename = f"pmns_p{self.params['p'].nbits()}_k{self.params['k']}_E{self.params['type']}.pkl"
        full_path = path / filename
        with open(full_path, 'wb') as f:
            pickle.dump(self, f)
            
            
            
            
    @classmethod
    def load(cls, pbits:int, k:int, Etype:int):
        path = OUTPUT_DIR_SAVES / f"pmns_p{pbits}_k{k}_E{Etype}.pkl"
        if not path.exists():
            raise FileNotFoundError(f"Path to file {path} doesn't seems to exist.")
        with open(path, 'rb') as f:
            return pickle.load(f)




    def add(self, name: str, element):
        if name in self.params:
            raise KeyError(f"Parameter '{name}' elready exist.")
        self.params[name] = self._ensure_list(element)
        
        
        
        
    def get(self, name: str, to_format: bool = False):
        val = self.params.get(name)
        
        if val is None:
            raise KeyError(f"Unknow parameter '{name}'.")

        if not to_format:
            return val

        n_limbs = self.params.get('n_limbs')
        
        is_list = isinstance(val, (list, tuple))
        is_nested = is_list and len(val) > 0 and isinstance(val[0], (list, tuple))
        is_tensor = is_nested and len(val[0]) > 0 and isinstance(val[0][0], (list, tuple))
        
        mpn_names = ['p', 'gamma', 'T_mat']
        if name in mpn_names:
            if is_nested: 
                return format.format_matrix_to_mpn(val, n_limbs)
            if hasattr(val, '_vector_'):
                return format.format_extension_field_to_mpn(val, n_limbs)
            return format.format_element_to_mpn(int(val), n_limbs)

        if is_tensor:
            return format.format_tensor_to_int64(val)
        if is_nested:
            return format.format_matrix_to_int64(val)
        if is_list:
            return format.format_list_to_int64(val)
            
        try:
            return format.format_to_int64(int(val))
        except (TypeError, ValueError):
            return str(val)
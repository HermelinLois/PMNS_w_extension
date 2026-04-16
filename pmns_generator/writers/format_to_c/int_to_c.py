# ==================================
# int_to_c.py
# Functions to write element composed
# with int64 or int128 in C
# ==================================

def cast_to_c(x:int) -> str:
    suffix = 'LL' if x < 0 else 'ULL'
    return str(x) + suffix


def cast_to_c128(x:int, fname:str) -> str:
    hi = x >> 64
    lo = x & (2**64 - 1)
    return f"{fname}({hi}{'LL' if x<0 else 'ULL'}, {lo}ULL)"


def format_list(l:list) -> str:
    return '{' + ', '.join(map(cast_to_c, l)) + '}'


def format_list128(l:list, fname:str) -> str:
    return '{' + ', '.join(map(lambda x : cast_to_c128(x, fname), l)) + '}'


def format_matrix(m, space:str =" "*4) -> str:
    lines = []
    for row in m:
        row_str = f"{space}{{{', '.join(str(x) + ('LL' if x < 0 else 'ULL') for x in row)}}}"
        lines.append(row_str)
    return '{\n' + ',\n'.join(lines) + '\n};'


def format_matrix128(m, fname:str, space:str =" "*4) -> str:
    mask = (1 << 64) - 1
    lines = []
    for row in m:
        formatted_elements = []
        for x in row:
            hi = x >> 64
            lo = x & mask
            h_suf = 'LL' if hi < 0 else 'ULL'
            formatted_elements.append(f"{fname}({hi}{h_suf}, {lo}ULL)")
        
        lines.append(f"{space}{{{', '.join(formatted_elements)}}}")
    return '{\n' + ',\n'.join(lines) + '\n};'


def format_element_to_BigInt(x, nb_chunks, phi_pow=64):
    decompose_x = []
    mask = (1 << phi_pow) - 1
    for chunk in range(nb_chunks):
        part = x & mask
        decompose_x.append(f"{part}ULL")
        x >>= phi_pow
    return  '{{' + ", ".join(decompose_x) + '}}'

def format_matrix_BigInt(m, nb_chunks, space:str =" "*4, phi_pow=64) -> str:
    lines = []
    for row in m:
        formatted_elements = [format_element_to_BigInt(x, nb_chunks, phi_pow) for x in row]
        row_str = space + "{\n" + space*2 + (",\n" + space*2).join(formatted_elements) + "\n" + space + "}"
        lines.append(row_str)

    return ",\n".join(lines).lstrip()

def format_element_to_mpn(x, nb_chunks, phi_pow=64):
    decompose_x = []
    mask = (1 << phi_pow) - 1
    for chunk in range(nb_chunks):
        part = x & mask
        decompose_x.append(f"{part}ULL")
        x >>= phi_pow
    return  '{' + ", ".join(decompose_x) + '}'

def format_matrix_to_mpn(m, nb_chunks, space:str =" "*4, phi_pow=64) -> str:
    lines = []
    for row in m:
        formatted_elements = [format_element_to_mpn(x, nb_chunks, phi_pow) for x in row]
        row_str = space + "{\n" + space*2 + (",\n" + space*2).join(formatted_elements) + "\n" + space + "}"
        lines.append(row_str)

    return ",\n".join(lines).lstrip()
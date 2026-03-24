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
    fmatrix = []
    for row in m:
        fmatrix.append(space + format_list(row))
    return '{\n' + ',\n'.join(fmatrix) + '\n};'

def format_matrix128(m,fname:str, space:str =" "*4) -> str:
    fmatrix = []
    for row in m:
        fmatrix.append(space + format_list128(row, fname))
    return '{\n' + ',\n'.join(fmatrix) + '\n};'
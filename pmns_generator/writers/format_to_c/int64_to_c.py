def cast_to_c(x:int) -> str:
    suffix = 'LL' if x < 0 else 'ULL'
    return str(x) + suffix

def format_list(l:list) -> str:
    return '{' + ', '.join(map(cast_to_c, l)) + '}'
    
def format_matrix(m, space:str =" "*4) -> str:
    fmatrix = []
    for row in m:
        fmatrix.append(space + format_list(row))
    return '{\n' + ',\n'.join(fmatrix) + '\n};'

def format_to_int64(x:int) -> str:
    suffix = 'LL' if x < 0 else 'ULL'
    return str(x) + suffix


def format_list_to_int64(l:list) -> str:
    return '{' + ', '.join(map(format_to_int64, l)) + '}'

def format_matrix_to_int64(m, space:str =" "*4) -> str:
    lines = []
    for row in m:
        row_str = f"{space}{{{', '.join(str(x) + ('LL' if x < 0 else 'ULL') for x in row)}}}"
        lines.append(row_str)
    return '{\n' + ',\n'.join(lines) + '\n}'


def format_element_to_mpn(x, nb_limbs, phi_pow=64) -> str:
    decompose_x = []
    mask = (1 << phi_pow) - 1
    for _ in range(nb_limbs):
        part = x & mask
        decompose_x.append(f"{part}ULL")
        x >>= phi_pow
    body = ", ".join(decompose_x)
    return "{" + body + "}"

def format_extension_field_to_mpn(x, nb_limbs, space = " "*4, phi_pow=64):
    decompose_x = [space + format_element_to_mpn(int(c), nb_limbs, phi_pow) for c in x._vector_()]
    return "{\n" + ",\n".join(decompose_x) + "\n}"
    
def format_matrix_to_mpn(m, nb_limbs, space:str =" "*4, phi_pow=64) -> str:
    lines = []
    for row in m:
        formatted_elements = [format_element_to_mpn(int(x), nb_limbs, phi_pow) for x in row]
        row_str = space + "{\n" + space*2 + (",\n" + space*2).join(formatted_elements) + "\n" + space + "}"
        lines.append(row_str)

    return "{\n" + ",\n".join(lines) + "\n}"


def format_tensor_to_mpn(t, nb_limbs, space=" " * 4, phi_pow=64):
    lines = []
    for mat in t:
        formatted_mat = format_matrix_to_mpn(mat, nb_limbs, space * 2, phi_pow)
        indented_mat = "\n".join((space * 2 + line) if line else line for line in formatted_mat.splitlines())
        row_str = space + "{\n" + indented_mat + "\n" + space + "}"
        lines.append(row_str)

    return "{\n" + ",\n".join(lines) + "\n}"


def format_tensor_to_int64(t, space=" " * 4):
    lines = []
    for mat in t:
        formatted_mat = format_matrix_to_int64(mat, space)
        lines.append(formatted_mat)

    return "{\n" + ",\n".join(lines) + "\n}"
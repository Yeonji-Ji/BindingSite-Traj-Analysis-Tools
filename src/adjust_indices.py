def adjust_index(file_path, indices):

    if file_path.endswith('.prmtop') or file_path.endswith('.nc'):
        indices = [idx - 1 for idx in indices]

    elif file_path.endswith('.pdb'):
        indices = [idx - 1 for idx in indices]

    else:
        print("Unsupported file format")

    return indices
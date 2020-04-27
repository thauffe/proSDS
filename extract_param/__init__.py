# get base frequencies from file
# input: path to file, substitution model
# output: list of base frequencies
def get_freq(sub_mod_dir, sub_mod):
    data = []
    data_start = False
    with open(sub_mod_dir, 'r') as file:
        for line in file:
            if data_start:
                if line.find('freq') != -1:
                    data.append(line.replace(' ', '').replace('\n', '').split('=')[-1])
                if line.find('Model') != -1:
                    break
            if line.find('Model = ' + sub_mod + '\n') != -1:
                data_start = True

    diff = 0
    for num in data:
        diff += float(num)
    data[0] = float(data[0]) + 1 - diff
    return data


# get substitution rate from file
# input: path to file, substitution model
# output: list of substitution rates
def get_R(sub_mod_dir, sub_mod):
    data = []
    for n in range(6):
        data.append("")
    data_start = False
    with open(sub_mod_dir, 'r') as file:
        for line in file:
            if data_start:

                if line.find('[CT]') != -1:
                    data[4] = line.replace(' ', '').replace('\n', '').split('=')[-1]
                if line.find('[AT]') != -1:
                    data[2] = line.replace(' ', '').replace('\n', '').split('=')[-1]
                if line.find('[GT]') != -1:
                    data[5] = line.replace(' ', '').replace('\n', '').split('=')[-1]
                if line.find('[AC]') != -1:
                    data[0] = line.replace(' ', '').replace('\n', '').split('=')[-1]
                if line.find('[CG]') != -1:
                    data[3] = line.replace(' ', '').replace('\n', '').split('=')[-1]
                if line.find('[AG]') != -1:
                    data[1] = line.replace(' ', '').replace('\n', '').split('=')[-1]
                if line.find('Model') != -1:
                    break

            if line.find('Model = ' + sub_mod + '\n') != -1:
                data_start = True

    return data


# get proportion of invariables and gamma shape from file
# input: path to file, substitution model
# output: list of invariables and gamma shape
def get_rateVarModel(sub_mod_dir, sub_mod):
    data = [0, 0]
    data_start = False

    with open(sub_mod_dir, 'r') as file:
        for line in file:
            if data_start:
                if line.find('p-inv') != -1:
                    data[0] = line.replace(' ', '').replace('\n', '').split('=')[-1]
                if line.find('gamma') != -1:
                    data[1] = line.replace(' ', '').replace('\n', '').split('=')[-1]
                if line.find('Model') != -1:
                    break
            if line.find('Model = ' + sub_mod + '\n') != -1:
                data_start = True

    return data


# get transition rates from file
# input: path to file, substitution model
# output: list transition rates
def get_ti_tv(sub_mod_dir, sub_mod):
    data = []
    data_start = False
    with open(sub_mod_dir, 'r') as file:
        for line in file:
            if line.find('Model = ' + sub_mod + '\n') != -1:
                data_start = True
            if data_start:
                if line.find('ti/tv') != -1:
                    data.append(line.replace(' ', '').replace('\n', '').replace(')', '').split('=')[-1])
                    break
    data.append(data[0])
    return data


# extracted params
# input:  path to file with params for the substitution model, birth rate, death rate, number of species
# output: list of params extracted from file
def extract_param(sub_mod_dir, sub_mod, birth_rate, death_rate, n_species):

    params = [[0], [0], [0]]

    sub_mod_mode = sub_mod.split('+')

    if sub_mod_mode[0] == 'GTR':
        params[1] = (get_R(sub_mod_dir, sub_mod))
        params[0] = (get_freq(sub_mod_dir, sub_mod))

    if sub_mod_mode[0] == 'HKY':
        params[1] = (get_ti_tv(sub_mod_dir, sub_mod))
        params[0] = (get_freq(sub_mod_dir, sub_mod))
    if len(sub_mod_mode) > 1:
        params[2] = (get_rateVarModel(sub_mod_dir, sub_mod))


    params.extend([birth_rate, death_rate, n_species])
    return params

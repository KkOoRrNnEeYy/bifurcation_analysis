import csv


def make_fieldnames(n_eq):
    fieldnames = [f'p{i}' for i in range(len(n_eq))]
    fieldnames.append('det')
    fieldnames.extend([f'sv{i}' for i in range(n_eq)])
    fieldnames.append('type')
    return fieldnames
    
def make_info(fieldnames, params, det, sv, pt):
    data = [*params, det, *sv, pt]
    info = {}
    for i in range(len(fieldnames)):
        info[fieldnames[i]] = data[i]
    return info

def create_csv(file, fieldnames):
    with open(file, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        
def write_csv(file, info):
    fieldnames = info.keys()
    with open(file, 'a') as f:
        csv_writer = csv.DictWriter(f, fieldnames=fieldnames)
        csv_writer.writerow(info)


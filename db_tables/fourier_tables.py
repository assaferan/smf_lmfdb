import sys
sys.path.append("/home/jean/code/lmfdb")
from lmfdb import db
from common_create_table import generate_table

def smf_qexp_short_col_type():
    cols = {}
    cols['smf_label'] = 'text'
    cols['nmax'] = 'smallint'
    cols['n1'] = 'smallint'
    cols['n2'] = 'smallint'
    cols['n12'] = 'smallint'
    cols['index'] = 'smallint'
    cols['coeff'] = 'integer[]'
    return cols

def smf_qexp_short_col_desc():
    desc = {}
    desc['smf_label'] = 'Label of the Siegel modular form'
    desc['nmax'] = 'Degree in q1, q2'
    desc['n1'] = 'Degree in q1'
    desc['n2'] = 'Degree in q2'
    desc['n12'] = 'Degree in q12'
    desc['index'] = 'Index of coefficient in case of vector-valued forms'
    desc['coeff'] = 'Coefficient encoded as a list of integers (coordinates in an integral basis of the Hecke field)'
    return desc

def smf_qexp_reps_col_type():
    cols = {}
    cols['smf_label'] = 'text'
    cols['qf'] = 'integer[]'
    cols['disc'] = 'integer'
    cols['index'] = 'smallint'
    cols['coeff'] = 'integer[]'
    return cols

def smf_qexp_reps_col_desc():
    desc = {}
    desc['smf_label'] = 'Label of the Siegel modular form'
    desc['qf'] = 'Reduced quadratic form encoded as a triple of integers'
    desc['disc'] = 'Discriminant of the quadratic form'
    desc['index'] = 'Index of coefficient in case of vector-valued forms'
    desc['coeff'] = 'Coefficient encoded as a list of integers (coordinates in an integral basis of the Hecke field)'
    return desc

def smf_qexp_short_header():
    header = "smf_label:nmax:n1:n2:n12:index:coeff\ntext:smallint:smallint:smallint:smallint:smallint:integer[]\n\n"
    return header

def smf_qexp_reps_header():
    header = "smf_label:qf:disc:index:coeff\ntext:integer[]:integer:smallint:integer[]\n\n"
    return header

#generate_table('smf_qexp_short', "Short, fully expanded q-expansions of Siegel modular forms", smf_qexp_short_col_type, smf_qexp_short_col_desc, label_col=None)

#generate_table('smf_qexp_reps', "q-expansions of Siegel modular forms indexed by reduced quadratic forms", smf_qexp_reps_col_type, smf_qexp_reps_col_desc, label_col=None)

def load_smf_qexp_reps():
    table = db.smf_qexp_reps
    aux_fname = "smf_qexp_reps_table.dat"
    header = smf_qexp_reps_header()
    with open(aux_fname, "w") as f:
        f.write(header)
    print_E4_qexp_reps(aux_fname)
    table.reload(aux_fname)
    return

def load_smf_qexp_short():
    table = db.smf_qexp_reps
    aux_fname = "smf_lmfdb/db_tables/smf_qexp_reps_table.dat"
    header = smf_qexp_short_header()
    with open(aux_fname, "w") as f:
        f.write(header)
    print_E4_qexp_short(aux_fname)
    table.reload(aux_fname)
    return

def E4_coefficients():
    coefs = {
        0: {
            (0, 0, 0): 1,
            (0, 0, 1): 240,
            (0, 0, 10): 272160,
            (0, 0, 11): 319680,
            (0, 0, 12): 490560,
            (0, 0, 13): 527520,
            (0, 0, 14): 743040,
            (0, 0, 15): 846720,
            (0, 0, 16): 1123440,
            (0, 0, 17): 1179360,
            (0, 0, 18): 1635120,
            (0, 0, 19): 1646400,
            (0, 0, 2): 2160,
            (0, 0, 20): 2207520,
            (0, 0, 21): 2311680,
            (0, 0, 22): 2877120,
            (0, 0, 23): 2920320,
            (0, 0, 24): 3931200,
            (0, 0, 25): 3780240,
            (0, 0, 3): 6720,
            (0, 0, 4): 17520,
            (0, 0, 5): 30240,
            (0, 0, 6): 60480,
            (0, 0, 7): 82560,
            (0, 0, 8): 140400,
            (0, 0, 9): 181680
        },
        3: {
            (1, 1, 1): 13440
        },
        4: {
            (1, 0, 1): 30240
        },
        7: {
            (1, 1, 2): 138240
        },
        8: {
            (1, 0, 2): 181440
        },
        11: {
            (1, 1, 3): 362880
        },
        12: {
            (1, 0, 3): 497280,
            (2, 2, 2): 604800
        },
        15: {
            (1, 1, 4): 967680,
            (2, 1, 2): 967680
        },
        16: {
            (1, 0, 4): 997920,
            (2, 0, 2): 1239840
        },
        19: {
            (1, 1, 5): 1330560
        },
        20: {
            (1, 0, 5): 1814400,
            (2, 2, 3): 1814400
        },
        23: {
            (1, 1, 6): 2903040,
            (2, 1, 3): 2903040
        },
        24: {
            (1, 0, 6): 2782080,
            (2, 0, 3): 2782080
        },
        27: {
            (1, 1, 7): 3279360,
            (3, 3, 3): 3642240
        },
        28: {
            (1, 0, 7): 4008960,
            (2, 2, 4): 5114880
        },
        31: {
            (1, 1, 8): 5806080,
            (2, 1, 4): 5806080
        },
        32: {
            (1, 0, 8): 5987520,
            (2, 0, 4): 7439040,
            (3, 2, 3): 5987520
        },
        35: {
            (1, 1, 9): 6531840,
            (3, 1, 3): 6531840
        },
        36: {
            (1, 0, 9): 7650720,
            (2, 2, 5): 7650720,
            (3, 0, 3): 8467200
        },
        39: {
            (1, 1, 10): 10644480,
            (2, 1, 5): 10644480,
            (3, 3, 4): 10644480
        },
        40: {
            (1, 0, 10): 9555840,
            (2, 0, 5): 9555840
        },
        43: {
            (1, 1, 11): 10039680
        },
        44: {
            (1, 0, 11): 13426560,
            (2, 2, 6): 16329600,
            (3, 2, 4): 13426560
        },
        47: {
            (1, 1, 12): 17418240,
            (2, 1, 6): 17418240,
            (3, 1, 4): 17418240
        },
        48: {
            (1, 0, 12): 15980160,
            (2, 0, 6): 19958400,
            (3, 0, 4): 15980160,
            (4, 4, 4): 20818560
        },
        51: {
            (1, 1, 13): 16208640,
            (3, 3, 5): 16208640
        },
        52: {
            (1, 0, 13): 18264960,
            (2, 2, 7): 18264960
        },
        55: {
            (1, 1, 14): 24192000,
            (2, 1, 7): 24192000,
            (4, 3, 4): 24192000
        },
        56: {
            (1, 0, 14): 23950080,
            (2, 0, 7): 23950080,
            (3, 2, 5): 23950080
        },
        59: {
            (1, 1, 15): 24312960,
            (3, 1, 5): 24312960
        },
        60: {
            (1, 0, 15): 28062720,
            (2, 2, 8): 35804160,
            (3, 0, 5): 28062720,
            (4, 2, 4): 35804160
        },
        63: {
            (1, 1, 16): 34974720,
            (2, 1, 8): 34974720,
            (3, 3, 6): 38707200,
            (4, 1, 4): 34974720
        },
        64: {
            (1, 0, 16): 31963680,
            (2, 0, 8): 39947040,
            (4, 0, 4): 41882400,
            (4, 4, 5): 31963680
        },
        67: {
            (1, 1, 17): 30360960
        },
        68: {
            (1, 0, 17): 38465280,
            (2, 2, 9): 38465280,
            (3, 2, 6): 38465280
        },
        71: {
            (1, 1, 18): 49351680,
            (2, 1, 9): 49351680,
            (3, 1, 6): 49351680,
            (4, 3, 5): 49351680
        },
        72: {
            (1, 0, 18): 42638400,
            (2, 0, 9): 42638400,
            (3, 0, 6): 47537280
        },
        75: {
            (1, 1, 19): 42349440,
            (3, 3, 7): 42349440,
            (5, 5, 5): 44029440
        },
        76: {
            (1, 0, 19): 49230720,
            (2, 2, 10): 59875200,
            (4, 2, 5): 49230720
        },
        79: {
            (1, 1, 20): 59996160,
            (2, 1, 10): 59996160,
            (4, 1, 5): 59996160
        },
        80: {
            (1, 0, 20): 59875200,
            (2, 0, 10): 74390400,
            (3, 2, 7): 59875200,
            (4, 0, 5): 59875200,
            (4, 4, 6): 74390400
        },
        83: {
            (1, 1, 21): 56246400,
            (3, 1, 7): 56246400
        },
        84: {
            (1, 0, 21): 63624960,
            (2, 2, 11): 63624960,
            (3, 0, 7): 63624960,
            (5, 4, 5): 63624960
        },
        87: {
            (1, 1, 22): 78382080,
            (2, 1, 11): 78382080,
            (3, 3, 8): 78382080,
            (4, 3, 6): 78382080
        },
        88: {
            (1, 0, 22): 67616640,
            (2, 0, 11): 67616640
        },
        91: {
            (1, 1, 23): 66528000,
            (5, 3, 5): 66528000
        },
        92: {
            (1, 0, 23): 84188160,
            (2, 2, 12): 107412480,
            (3, 2, 8): 84188160,
            (4, 2, 6): 107412480
        },
        95: {
            (1, 1, 24): 101606400,
            (2, 1, 12): 101606400,
            (3, 1, 8): 101606400,
            (4, 1, 6): 101606400,
            (5, 5, 6): 101606400
        },
        96: {
            (1, 0, 24): 91808640,
            (2, 0, 12): 114065280,
            (3, 0, 8): 91808640,
            (4, 0, 6): 114065280,
            (4, 4, 7): 91808640,
            (5, 2, 5): 91808640
        },
        99: {
            (1, 1, 25): 85276800,
            (3, 3, 9): 95074560,
            (5, 1, 5): 85276800
        },
        100: {
            (1, 0, 25): 93774240,
            (2, 2, 13): 93774240,
            (5, 0, 5): 97554240
        }
    }
    return coefs

def print_E4_qexp_reps(filename):
    smf_label = '2.K.1.4.0.a.a'
    coefs = E4_coefficients()
    with open(filename, "a") as f:
        for disc in coefs.keys():
            for tup in (coefs[disc]).keys():
                coef = coefs[disc][tup]
                line = "{}:[{},{},{}]:{}:{}:[{}]".format(smf_label, tup[0], tup[1], tup[2], disc, 0, coef)
                line = line.replace("[", "{")
                line = line.replace("]", "}")
                f.write(line)
                f.write("\n")
    return

def print_E4_qexp_short(filename):
    smf_label = '2.K.1.4.0.a.a'
    with open(filename, "a") as f:
        #label,nmax,n1,n2,n12,coeff
        f.write("{}:0:0:0:0:{{1}}\n".format(smf_label))
        f.write("{}:1:1:0:0:{{240}}\n".format(smf_label))
        f.write("{}:1:0:0:1:{{240}}\n".format(smf_label))
        f.write("{}:1:1:2:1:{{240}}\n".format(smf_label))
        f.write("{}:1:1:1:1:{{13440}}\n".format(smf_label))
        f.write("{}:1:1:0:1:{{30240}}\n".format(smf_label))
        f.write("{}:1:1:-1:1:{{13440}}\n".format(smf_label))
        f.write("{}:1:1:-2:1:{{240}}\n".format(smf_label))


import signal
from lmfdb import db

# populating from samples
# start with only one

sampls = db.smf_samples.search()
sam = [s for s in sampls]

# at the moment these are the only families we have in the samples here
def get_newform_family(collection):
    if collection == 'Kp':
       return 'K'
    else:
       return 'F'

def get_newform_level(collection, name):
    if collection == 'Kp':
       return eval(name.split('_')[2])
    else:
       return 1

def get_newform_weight(collection, weight, degree):
    if collection == 'Sp4Z_2':
       return [2 + weight, weight]
    elif collection == 'Sp4Z_j':
       raise RuntimeError("Don't know what j is!")
    else:
       return [weight for i in [1..degree]]

def dark_magic(z,g):
    M = matrix([list(z[g+i*(i-1)/2:g+i*(i+1)/2]) + [z[i]/2] + [0 for j in range(g-1-i)] for i in range(g)]) 
    A = M + M.transpose() 
    return A 

# this doesn't really work
def handler(signum, frame):
   raise Exception("time's up!")

signal.signal(signal.SIGALRM, handler)

def create_dict_from_sample(s, entries):
    e = dict()                                                                
    e['id'] = s['id_link'] 
    e['family'] = get_newform_family(s['collection'][0])
    e['degree'] = s['degree']                                                    
    e['level'] = get_newform_level(s['collection'][0], s['name'])
    e['dim'] = s['fdeg']                                                          
    e['weight'] = get_newform_weight(s['collection'][0], s['weight'], s['degree'])
    e['level_is_prime'] = is_prime(e['level'])                                
    e['level_is_prime_power'] = is_prime_power(e['level'])                    
    e['level_is_square'] = is_square(e['level'])                              
    e['level_is_squarefree'] = is_squarefree(e['level'])                      
    e['level_primes'] = prime_divisors(e['level'])
    e['weight_alt'] = [e['weight'][-1]] + [e['weight'][-i-1] - e['weight'][-i] for i in range(1, e['degree'])]
    last_nz = max([-1]+[i for i in range(e['degree']) if e['weight_alt'][i] != 0])
    e['weight_alt'] = [e['weight_alt'][i] for i in range(last_nz+1)]
    e['char_orbit_label'] = 'a'
    e['space_label'] = '.'.join([str(e['degree']), e['family'], str(e['level'])] + [str(l) for l in e['weight_alt']] + [e['char_orbit_label']])
    others = [f for f in entries if f['space_label'] == e['space_label']]
    e['hecke_orbit'] = 1 + len(others)
    e['label'] = e['space_label'] + '.' + chr(ord('a') + (e['hecke_orbit']-1))
    e['char_orbit_index'] = 1
    e['hecke_orbit_code'] = e['degree'] + (ord(e['family'])<<8) + (e['level']<<12)+(e['weight_alt'][0] << 20) + ((e['char_orbit_index']-1)<<36) + ((e['hecke_orbit']-1)<<52)
    e['is_polredabs'] = false
    e['char_degree'] = 1                                                      
    e['char_is_minimal'] = true                                               
    e['char_conductor'] = 1                                                   
    e['prim_orbit_index'] = 1                                                 
    e['char_is_real'] = true
    e['char_parity'] = 1
    _.<x> = PolynomialRing(RationalField())
    F = sage_eval(s['field'], locals={'x':x})
    def trace(F, F_elt):
        if F.degree() == 1:
           return F_elt
        else:
           return F(F_elt).trace()
    a = F.0
    pol = F.defining_polynomial()
    e['field_poly'] = pol.coefficients(sparse=False)
#    signal.alarm(60)
#    try:
#      e['field_disc'] = F.discriminant()
#    except exc:
#      e['field_disc'] = ""
#    signal.alarm(0)
#    if e['id'] in [57, 66, 67, 70, 71, 73, 78, 85, 86, 87, 88, 101, 102, 104, 108, 109, 127, 128, 129]:
    if F.degree() > 5 and e['id'] <= 129:
      e['field_disc'] = ""
    else:
      e['field_disc'] = F.discriminant()
    e['nf_label'] = "NULL"
    if e['field_disc'] != "":
      fields = [t for t in db.nf_fields.search({'disc_abs' : abs(e['field_disc'])})]
      if len(fields) > 0:
         field = fields[0]
         e['nf_label'] = field['label']
    e['char_order'] = 1
    e['related_objects'] = []
    e['embedded_related_objects'] = []
    e['field_poly_is_cyclotomic'] = False
    e['field_poly_is_real_cyclotomic'] = False
#    if F.degree() == 1 or ((e['id'] not in [70,85,86,87,101,108]) and F.is_abelian()):
    if (F.degree() == 1) or ((F.degree() <= 5) and (F.is_abelian())):
       if F.degree() == 1:
          cond = 1
       else:
          cond = F.conductor()
       e['field_poly_is_cyclotomic'] = (F.degree() == euler_phi(cond)) and (pol == cyclotomic_polynomial(cond))
       if (2*F.degree() == euler_phi(cond)):
          K.<zeta> = CyclotomicField(cond)
          omega = zeta + zeta^(-1)
          e['field_poly_is_real_cyclotomic'] = (omega.min_poly() == pol)
    e['analytic_rank'] = -1
    e['field_poly_root_of_unity'] = 0
    if e['field_poly_is_cyclotomic'] or e['field_poly_is_real_cyclotomic']:
       e['field_poly_root_of_unity'] = cond
    e['level_radical'] = product(e['level_primes'])
    e['analytic_rank_proved'] = False
    if e['field_disc'] == "" or e['id'] == 74:
      e['field_disc_factorization'] = ""
    else:
      e['field_disc_factorization'] = list(e['field_disc'].factor())
      e['field_disc_factorization'] = str(e['field_disc_factorization']).replace('(', '{').replace(')', '}')
    e['relative_dim'] = e['dim']
    fc_data = [x['data'] for x in db.smf_fc.search({'owner_id' : str(e['id'])})]
    e['qexp_display'] = ('+'.join(['+'.join([x[1] + 'q^{' + latex(dark_magic(eval(x[0]),e['degree'])) + '}' for x in t.items()]) for t in fc_data[:10]])).replace('+-','-') + '+\\cdots'
    ev_data = [x for x in db.smf_ev.search({'owner_id' : e['id']})]
    tr = [[d['data'] for d in ev_data if d['index'] == p] for p in [2,3,5,7]]
    e['trace_display'] = [trace(F, sage_eval(t[0], locals={F.variable_name():F.0})) for t in tr if len(t) > 0]
    # we only display integral traces, which might not happen for Eisenstein forms
    good_indices = []
    for i in range(len(e['trace_display'])):
       x = e['trace_display'][i]
       if (x.is_integral()):
          good_indices.append(i)
    e['trace_display'] = [e['trace_display'][i] for i in good_indices]
    max_idx = max([1] + [d['index'] for d in ev_data])
    tr_all = [[d['data'] for d in ev_data if d['index'] == n] for n in range(2,max_idx+1)]
    e['traces'] = [1] + [trace(F, sage_eval(t[0], locals={F.variable_name():F.0})) for t in tr_all if len(t) > 0]
    for i in range(len(e['traces'])):
       x = e['traces'][i]
       if (x.is_integral()):
       	  good_indices.append(i)
    e['traces'] = [e['traces'][i] for i in good_indices]
    e['trace_hash'] = hash(str(e['traces'])) % 2^63
    e['qexp_display'] = e['qexp_display'].replace('\n','').replace("\\", "\\\\")
    e['is_cuspidal'] = False
    if 'cusp form' in s['type']:
       e['is_cuspidal'] = true
    lift = [w.split()[0] for w in s['type'].split(',') if 'lift' in w]
    if len(lift) > 0:
       e['lift_type'] = lift[0]
    else:
       e['lift_type'] = "None"
    return e

def write_data(entries):
    column_names = '|'.join(db.smf_newforms.col_type.keys())
    column_types = '|'.join(db.smf_newforms.col_type.values())
    e_data = []
    for e in entries:
        e_datum = '|'.join([str(e[k]) for k in db.smf_newforms.col_type.keys()])
        e_datum = e_datum.replace('[', '{').replace(']','}')
        e_data.append(e_datum)
    write_data = "\n".join([column_names, column_types, ""] + e_data)
    f = open("../smf/table_data.dat", "w")
    f.write(write_data)
    f.close()
    return

entries = []
for s in sam:
  print("s['id_link'] = ", s['id_link'])
  e = create_dict_from_sample(s, entries)
  entries.append(e)
write_data(entries)                                                                  
db.smf_newforms.reload("../smf/table_data.dat", null="")
db.smf_newforms.cleanup_from_reload()
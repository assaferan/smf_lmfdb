db.create_table("smf_newspaces", {"text" : ["label", "family"], "smallint": "degree", "integer" : ["level", "dim"]}, "label", table_description = "Siegel newspaces", col_description = {"label" : "Label g.C.N.w.a of this newspace", "family" : "Family of congruence subgroups it belongs to", "degree" : "Degree g of the Siegel modular group", "level" : "Level of the arithmetic subgroup", "dim" : "Dimension of the space"})

db.smf_newspaces.column_description("label", "Label g.C.N.w.a of this newspace")

db.smf_newspaces.column_description("family", "Family of arithmetic subgroups ('F' = full, 'K' = paramodular, 'S' = Siegel, 'C' - principal)")                         

db.smf_newspaces.column_description("degree", "Degree g of this newform (automorphic with repsect to Sp(2g, Q))")                                    

db.smf_newspaces.column_description("level", "Level N in the family")

db.smf_newspaces.add_column('weight', 'smallint[]', "Weight of this newform (highest weight of the corresponding irreducible representation of GL(g))" )                

db.smf_newspaces.add_column('level_is_prime', 'boolean', 'true if N is prime (1 is not prime)')                                                      

db.smf_newspaces.add_column('level_is_prime_power', 'boolean', 'true if N is a prime power (1 is not a prime power)')                                

db.smf_newspaces.add_column('level_is_square', 'boolean', 'true if N is square')

db.smf_newspaces.add_column('level_is_squarefree', 'boolean', 'true if N is squarefree')                                                             

db.smf_newspaces.add_column('level_primes', 'integer[]', 'sorted list of prime divisors of N')                                                       

db.smf_newspaces.column_description('dim', 'the Q-dimension of this newspace')
db.smf_newspaces.add_column('cusp_dim', 'integer', 'the Q-dimension of the space of cusp forms of this level, weight and character')
db.smf_newspaces.add_column('mf_dim', 'integer', 'the Q-dimension of the space of modular forms of this level, weight and character')
db.smf_newspaces.add_column('mf_new_dim', 'integer', 'the Q-dimension of the space of new modular forms of this level, weight and character')
db.smf_newspaces.add_column('eis_dim', 'integer', 'the Q-dimension of the space of Eisenstein forms of this level, weight and character')
db.smf_newspaces.add_column('eis_new_dim', 'integer', 'the Q-dimension of the space of new Eisenstein forms of this level, weight and character')

db.smf_newspaces.add_column('num_forms', 'smallint', 'number of Hecke orbits (each corresponds to a Galois conjugacy class of modular forms)')

db.smf_newspaces.add_column('char_order', 'integer', 'the order of the character chi')

db.smf_newspaces.add_column('level_radical', 'integer', 'product of prime divisors of N')                                                            

db.smf_newspaces.add_column('trace_display', 'numeric[]', 'list of integer traces tr(a_2), tr(a_3), tr(a_5), tr(a_7), only set when dim > 0 and not yet computed in every case.')

db.smf_newspaces.add_column('traces', 'numeric[]', 'integer coefficients a_n of the trace form (sum of all newforms in the space) for n from 1 to 1000, only set when dim > 0 and not yet computed in every case.')                         

db.smf_newspaces.add_column('char_parity', 'smallint', '1 for even, -1 for odd')    

db.smf_newspaces.add_column('char_degree', 'integer', 'Degree of the (cyclotomic) character field')                                                  

db.smf_newspaces.add_column('char_is_minimal', 'boolean', "true if the character chi is {{KNOWL('character.dirichlet.minimal','minimal')}}")         

db.smf_newspaces.add_column('char_conductor', 'integer', 'Conductor of the Dirichlet character chi of this newform')                                 

db.smf_newspaces.add_column('prim_orbit_index', 'smallint', 'char_orbit for the primitive version of this character')                                

db.smf_newspaces.add_column('char_is_real', 'boolean', 'true if the character takes only real values (trivial or quadratic)')                        

db.smf_newspaces.add_column('char_orbit_index', 'smallint', 'ordinal i identifying the Galois orbit of the character of this newform (base26 encoded in the newform label / character orbit label)')                           

db.smf_newspaces.add_column('char_orbit_label', 'text', 'base26-encoding of char_orbit_index-1')                                                     

db.smf_newspaces.add_column('weight_alt', 'integer[]', 'the tuple lambda_g, lambda_{g-1} - lambda_g, ..., omitting trailing zeros (when g = 2, these are (k,j) and for scalar valued we just get k) ')
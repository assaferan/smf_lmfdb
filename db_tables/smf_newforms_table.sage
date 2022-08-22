db.create_table("smf_newforms", {"text" : ["label", "family"], "smallint": "degree", "integer" : ["level", "weight", "dim"]}, "label", table_description = "Siegel newforms", col_description = {"label" : "Label C.N.k.j.a.x of this newform", "family" : "Family of congruence subgroups it belongs to", "degree" : "Degree g of the Siegel modular group", "level" : "Level of the congruence subgroup", "weight" : "Weight of the modular form", "dim" "Degree of the number field over which the form is defined"})

db.smf_newforms.column_description("label", "Label g.C.N.w.a.x of this newform")

db.smf_newforms.column_description("family", "Family of arithmetic subgroups ('F' = full, 'K' = paramodular, 'S' = Siegel, 'C' - principal)")                         

db.smf_newforms.column_description("degree", "Degree g of this newform (automorphic with repsect to Sp(2g, Q))")                                    

db.smf_newforms.column_description("level", "Level N in the family")

db.smf_newforms.drop_column('weight')

db.smf_newforms.add_column('weight', 'smallint[]', "Weight of this newform (highest weight of the corresponding irreducible representation of GL(g))" )                

db.smf_newforms.add_column('level_is_prime', 'boolean', 'true if N is prime (1 is not prime)')                                                      

db.smf_newforms.add_column('level_is_prime_power', 'boolean', 'true if N is a prime power (1 is not a prime power)')                                

db.smf_newforms.add_column('level_is_square', 'boolean', 'true if N is square')

db.smf_newforms.add_column('level_is_squarefree', 'boolean', 'true if N is squarefree')                                                             

db.smf_newforms.add_column('level_primes', 'integer[]', 'sorted list of prime divisors of N')                                                       

db.smf_newforms.add_column('space_label', 'text', 'label g.C.N.w.a of the newspace containing this newform')                                        

db.smf_newforms.add_column('hecke_orbit', 'integer', '(X) An integer that is encoded into x in the label via 1=a, 2=b, 26=z, 27=ba, 28=bb.  Note the shift: the letter is the Cremona code for X-1.')                         

db.smf_newforms.add_column('hecke_orbit_code', 'bigint', 'encoding of the tuple (g.C.N.w.i.x) into 64 bits, used in eigenvalue tables.  g + (ord(C)<<8) + (N<<12) + (w[0]<<20) + (w[1]<<28) + ((i-1)<<36) + ((X-1)<<52).')

db.smf_newforms.column_description('dim', 'the dimension of this newform')

db.smf_newforms.add_column('is_polredabs', 'boolean', "whether the polynomial has been reduced by Pari's `polredabs`")                              

db.smf_newforms.add_column('nf_label', 'text', 'LMFDB label for the corresponding number field (can be NULL)')                                      

db.smf_newforms.add_column('trace_hash', 'bigint', 'appropriate linear combination of the a_{p,i} between 2^12 and 2^13')                           

db.smf_newforms.add_column('analytic_rank', 'smallint', 'the analytic rank of the L-function of the embedded newforms in this newform orbit')       

db.smf_newforms.add_column('qexp_display', 'text', 'latexed string for display on search page results')                                             

db.smf_newforms.add_column('char_order', 'integer', 'the order of the character chi')                                                               

db.smf_newforms.add_column('related_objects', 'text[]', 'list of text URLs of related objects (e.g. elliptic curve isogeny class, Artin rep, ...), e.g. ["EllipticCurve/Q/11/a"]')                                            

db.smf_newforms.add_column('embedded_related_objects', 'text[]', 'list of lists of text URLs of related objects (e.g. Artin reps), indexed by embedding_m (so first entry is a list of friends for the first embeddeded newform)')                                                                      
db.smf_newforms.add_column('field_disc', 'numeric', 'discriminant of the coefficient field, if known')                                              

db.smf_newforms.add_column('field_poly_is_cyclotomic', 'boolean', "true if field_poly is a cylcotomic polynomial (the field might be Q(zeta_n) even when this flage is not set if we haven't chosen a cyclotomic polynomial to define it)")               

db.smf_newforms.add_column('field_poly', 'numeric[]', 'list of integers giving defining polynomial for the Hecke field (standard Sage order of coefficients)')         

db.smf_newforms.add_column('field_poly_is_real_cyclotomic', 'boolean', "true if field_poly is the minimal polynomial of zeta_n + zeta_n^-1 for some n (the field might be Q(zeta_n)^+ even when this flage is not set if we haven't chosen a cyclotomic polynomial to define it)")                      

db.smf_newforms.add_column('field_poly_root_of_unity', 'integer' , 'the value of n if either field_poly_is_cylotomic of field_poly_is_real_cyclotomic is set')         

db.smf_newforms.add_column('level_radical', 'integer', 'product of prime divisors of N')                                                            

db.smf_newforms.add_column('analytic_rank_proved', 'boolean', 'true if the analytic rank has been proved')                                          

db.smf_newforms.add_column('field_disc_factorization', 'numeric[]', 'factorization of field discriminant stored as ordered list of pairs [p,e]')    

db.smf_newforms.add_column('relative_dim', 'integer', 'the Q(chi)-dimension of this Hecke orbit (=dim/char_degree)')                                

db.smf_newforms.add_column('trace_display', 'numeric[]', 'list of the first four a_{p,1} traces for display on search page results')

db.smf_newforms.add_column('traces', 'numeric[]', 'full list of integer traces tr(a_{n,1}) for n from 1 to 1000 (or more)')                         

db.smf_newforms.add_column('char_parity', 'smallint', '1 for even, -1 for odd')    

db.smf_newforms.add_column('char_degree', 'integer', 'Degree of the (cyclotomic) character field')                                                  

db.smf_newforms.add_column('char_is_minimal', 'boolean', "true if the character chi is {{KNOWL('character.dirichlet.minimal','minimal')}}")         

db.smf_newforms.add_column('char_conductor', 'integer', 'Conductor of the Dirichlet character chi of this newform')                                 

db.smf_newforms.add_column('prim_orbit_index', 'smallint', 'char_orbit for the primitive version of this character')                                

db.smf_newforms.add_column('char_is_real', 'boolean', 'true if the character takes only real values (trivial or quadratic)')                        

db.smf_newforms.add_column('char_orbit_index', 'smallint', 'ordinal i identifying the Galois orbit of the character of this newform (base26 encoded in the newform label / character orbit label)')                           

db.smf_newforms.add_column('char_orbit_label', 'text', 'base26-encoding of char_orbit_index-1')                                                     

db.smf_newforms.add_column('char_values', 'jsonb', 'quadruple <N,n,u,v> where N is the level, n is the order of the character, u is a list of generators for the unit group of Z/NZ, and v is a corresponding list of integers for which chi(u[i]) = zeta_n^v[i]')

db.smf_newforms.add_column('weight_alt', 'integer[]', 'the tuple lambda_g, lambda_{g-1} - lambda_g, ..., omitting trailing zeros (when g = 2, these are (k,j) and for scalar valued we just get k) ')

db.smf_newforms.add_column('is_cuspidal', 'boolean', 'true if this is a cusp form')

db.smf_newforms.add_column('lift_type', 'text', 'If not a lift, states None')
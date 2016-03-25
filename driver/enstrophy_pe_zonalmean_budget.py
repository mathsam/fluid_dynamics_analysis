enstrophy_budget = EnstrophyTotal(psic)
ens_adv = enstrophy_budget.advection()
ens_gen = enstrophy_budget.enstrophy_generation()
ens_gen_ms = enstrophy_budget.enstrophy_generation_shear()
##
pe_budget = PETotal(psic, 500)
pe_adv = pe_budget.pe_adv()
bc2bt  = pe_budget.bc2bt()
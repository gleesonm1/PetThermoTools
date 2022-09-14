using MAGEMin_C

function satPhase(P_kbar, T_C, bulk)

		gv, DB = init_MAGEMin();
		sys_in = "wt";

		gv.verbose = -1;

		new_bulk = bulk/sum(bulk);
		
		out = point_wise_minimization(P_kbar, T_C, new_bulk, gv, DB, sys_in);
		Phase = out.ph
		Oxides = out.oxides
		if "liq" in Phase
			Liq_index = findfirst(x -> occursin("liq", x), out.ph)
			Liq_Comp = out.SS_vec[Liq_index].Comp_wt
			Liq_Frac = out.ph_frac_wt[Liq_index]
		
			finalize_MAGEMin(gv,DB);
			Ret = Dict("Phase" => Phase, "Oxides" => Oxides, "Liq_Comp" => Liq_Comp, "Liq_Frac" => Liq_Frac)
			return Ret
		else
			finalize_MAGEMin(gv, DB);
			Ret = Dict("Phase" => Phase, "Oxides" => Oxides)
			return Ret
		end
	end
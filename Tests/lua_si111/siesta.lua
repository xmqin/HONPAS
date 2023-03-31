
--[[
   siesta communication functions
--]]

is_broyden = false
is_linear = false
function siesta_comm()

   if siesta.state == siesta.SCF_LOOP then

      -- First retrieve the dDmax
      siesta.receive({'SCF.dD', 'SCF.dH'})
      siesta.SCF.dH = siesta.SCF.dH / siesta.Units.eV

      -- Check values
      if siesta.SCF.dD < 0.005 and siesta.SCF.dH < 0.005 and
        not is_linear then
	 
	 siesta.receive({'SCF.Mixer.Switch'})
	 siesta.SCF.Mixer.Switch = 'Linear'
	 siesta.send({'SCF.Mixer.Switch'})
	 is_linear = true
	 
      elseif siesta.SCF.dD < 0.02 and siesta.SCF.dH < 0.1 and
        not is_broyden then

	 siesta.receive({'SCF.Mixer.Switch'})
	 siesta.SCF.Mixer.Switch = 'Broyden'
	 siesta.send({'SCF.Mixer.Switch'})
	 is_broyden = true
      end
      
   end
   
end


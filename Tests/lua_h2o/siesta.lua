
--[[
   siesta communication functions
--]]

-- Global container for data
data = { tot_itt = 0, wmix = {} }

function siesta_comm()

   --[[
      Retrieve siesta information.

      In this case we only request the coordinates,
      the forces and the mixing weight
   --]]
   local get_tbl = {"SCF.Mixer.Weight","geom.xa","geom.fa"}
   local ret_tbl = {}

   -- Do the actual communication with fortran
   siesta.receive(get_tbl)

   ---[[
      -- You can even write to other files while running 
      local fh = assert(io.open("lua_out_".. siesta.Node ..".out","a"))
      fh:write("Current state: "..siesta.state.."\n")
      fh:write("Current wmix: "..siesta.SCF.Mixer.Weight.."\n")
      fh:close()
   --]]

   -- Print out what is available
   if siesta.state == siesta.INITIALIZE then
      -- Change variables if needed
      siesta:print("...After initialization...")
   end

   if siesta.state == siesta.INIT_MD then
      siesta:print("...Right before entering the SCF loop...")
      local fh = assert(io.open("lua_out_".. siesta.Node ..".out","a"))
      write_mat(fh,"xa",siesta.geom.xa)
      fh:close()
      ret_tbl = init_md(siesta)
   end

   if siesta.state == siesta.SCF_LOOP then
      siesta:print("...At start of SCF...")
      ret_tbl = scf(siesta)
   end

   if siesta.state == siesta.FORCES then
      siesta:print("...At atomic movement...")
      ret_tbl = siesta_move(siesta)
   end
   
   if siesta.state == siesta.ANALYSIS then
      siesta:print("...At analysis...")
      --[[
	 Uncomment the below function to make lua
	 print a list of available variables
	 
	 This is in this example called at the analysis
	 step
      --]]
      siesta.print_allowed()
   end

   --[[
      The ret_tbl returns only those
      variables back to siesta.

      In this way you can retrieve the coordinates
      and the forces, but only for post-processing
      and for creation special output formats.

      For MD calculations you may even retrieve
      forces and coordinates to _only_ update
      the coordinates, thus only sending back
      the coordinates.
      
      Needless to say, less communication, more
      speed, yet the overhead is minimal.
   --]]
   siesta.send(ret_tbl)
end

function init_md(siesta)
   
   -- One can check what the meshcutoff is etc.
   return {}
   
end

function scf(siesta)

   -- just create a reference (easier to read)
   SCF = siesta.SCF

   --[[
      Easy way to retrieve current SCF step
      siesta.SCF.Iteration also exists for this purpose
      siesta.SCF.Iteration gets reset after each MD move
      the below does not.
   --]]
   data.tot_itt = data.tot_itt + 1
   local i = data.tot_itt
   data.wmix[i] = SCF.Mixer.Weight
   
   if i > 1 and SCF.Mixer.Weight < 0.5 then
      --[[
	 We increase the mixing weight until a limit
	 of 0.5
      --]]
      SCF.Mixer.Weight = SCF.Mixer.Weight + 0.01
   end
   
   -- Only return mixing weight here
   return {"SCF.Mixer.Weight"}
end
   
function siesta_move(siesta)
   fa = siesta.geom.fa

   --[[
      This constraint is the same
      as 
        position from 2 to 3
      which shows the versatility of the Lua
      interface.
      
      Note that both Geometry.Constraints
      and the Lua calls work together,
      this way you can use Geometry.Constraints
      for the simpler constraints and Lua for
      more complex ones.
   --]]
   for i = 1 , #fa do
      if i > 1 then
	 fa[i][1] = 0.
	 fa[i][2] = 0.
	 fa[i][3] = 0.
      end
   end

   -- After MD step, we reset mixing-weight
   siesta.SCF.Mixer.Weight = data.wmix[1]

   -- Send back mixing weight and the corrected forces
   return {"SCF.Mixer.Weight","geom.fa"}
end

--[[
   Helper function to print out a matrix to
   a file-handle.
--]]
function write_mat(fh,name,mat)
   fh:write("Variable: "..name.."\n")
   for ia,xyz in pairs(mat) do
      if type(xyz) == "table" then
	 a = "  "
	 for _,x in pairs(xyz) do a = a .. " " .. x end
	 a = a .. "\n"
	 fh:write(a)
      end
   end
end


-- Load the flos module.
-- You have to ensure that flos is in the LUA_PATH
-- environment variable.
-- See:
--   https://github.com/siesta-project/flos
-- for details
local flos = require "flos"

relax = flos.LBFGS()

-- Define the communication function for
-- SIESTA
function siesta_comm()

   -- Perform relaxation and return information to
   -- SIESTA
   relax:SIESTA(siesta)
end

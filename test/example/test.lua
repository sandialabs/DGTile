-- this is an example application input
-- 

input = {

  simulation = {
    name = "advection",
    polynomial_order = 1,
    quadrature_rule = 2,
    tensor_product = true
  }

  mesh = {
    X = {num_blocks = 2, num_cells = 10, min = 0., max = 1.},
    Y = {num_blocks = 2, num_cells = 10, min = 0., max = 1.},
  },

  initial_conditions = {
  }

}

dgtile.run(input)

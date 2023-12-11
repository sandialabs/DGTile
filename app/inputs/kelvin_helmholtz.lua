density0 = function(x,y,z)
  if (math.abs(y) > 0.25) then
    return 1.
  else
    return 2.
  end
end

velocity0 = function(x,y,z)
  vx = 0.
  exp1 = math.exp(-(y-0.25)*(y-0.25) / 0.05)
  exp2 = math.exp(-(y+0.25)*(y+0.25) / 0.05)
  sin = math.sin(4.*math.pi*x)
  vy = 0.1*sin*(exp1 + exp2)
  if (math.abs(y) > 0.25) then
    vx = -0.5
  else
    vx = 0.5
  end
  return vx, vy, 0.
end

input = {

  name = 'kelvin_helmholtz',

  time = {
    cfl = 0.9,
    end_time = 1.,
    to_terminal = {{kind='at_step_periodically', frequency=1}},
    integrator = 'ssprk2'
  },

  basis = {
    polynomial_order = 1,
    quadrature_rule = 2,
    tensor_product = true
  },

  mesh = {
    X = {num_blocks = 1, num_cells = 80, min = -0.5, max = 0.5, periodic = true},
    Y = {num_blocks = 1, num_cells = 80, min = -0.5, max = 0.5, periodic = true}
  },

  hydro = {
    gamma = 1.4,
    density = density0,
    velocity = velocity0,
    pressure = 2.5
  },

  vtk = {
    when = {{kind='at_time_periodically', frequency=0.1}}
  }

}

dgtile.run(input)

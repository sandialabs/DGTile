-- build standard basis, where quadrature rule is set to p+1
for dim = 1,3,1 do
  for p = 0,3,1 do
    q = p+1
    basis = {
      dimension = dim,
      polynomial_order = p,
      quadrature_rule = q,
      tensor_product = true
    }
    dgt_check_basis(basis)
  end
end

-- build under-integrated linear basis
for dim = 1,3,1 do
  p = 1
  q = 1
  basis = {
    dimension = dim,
    polynomial_order = p,
    quadrature_rule = q,
    tensor_product = false
  }
  dgt_check_basis(basis)
end

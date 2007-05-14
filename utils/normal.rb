SQ2PI = Math.sqrt(2 * Math::PI)

def norm_dist(z)
  z = z.to_f
  z2 = z * z
  t = q = z * Math.exp(-0.5 * z2) / SQ2PI
  3.step(199, 2) do |i|
    prev = q
    t *= z2 / i
    q += t
    if q == prev
      return q
    end
  end
  z > 0 ? 0.5 : -0.5
end

def newton_a(y, ini, epsilon = 0.00000000001, limit = 30)
  x = ini
  limit.times do |i|
    prev = x
    f, df = yield(prev)
    x = (y - f)/df + prev
    break if (x - prev).abs < epsilon
  end
  x
end

def p_norm(y)
  newton_a(y, 0.0) {|x| [norm_dist(x), Math.exp(-0.5 * x**2) / SQ2PI]}
end



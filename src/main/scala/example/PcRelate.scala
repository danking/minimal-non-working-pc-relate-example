import breeze.linalg._

object PcRelate {
  def apply() {
    val g = new DenseMatrix(4,8,Array[Double](0,0,0,0, 0,0,1,0, 0,1,0,1, 0,1,1,1, 1,0,0,0, 1,0,1,0, 1,1,0,1, 1,1,1,1))
    val v = new DenseMatrix(4,2, Array[Double](0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0))

    val oneV = DenseMatrix.horzcat(new DenseMatrix(4,1, Array[Double](1,1,1,1)), v)

    val beta = (inv(oneV.t * oneV) * oneV.t * g).t

    val mu_si = (beta * oneV.t) / 2.0
    val mu_si_clipped = mu_si.map(x => if (x <= 0.0) Double.MinPositiveValue else (if (x >= 1.0) 1 - Double.MinPositiveValue else x))

    val g2mu = g.t - (2.0 * mu_si_clipped)
    val numer = g2mu.t * g2mu
    val variance = mu_si_clipped :* (1.0 - mu_si_clipped)
    val denom = 4.0 :* (variance.t.map(math.sqrt) * variance.map(math.sqrt))
    val phi_hat = numer :/ denom

    println("phi")
    println(phi_hat)
    println("")

    val mu_is_clipped = mu_si_clipped.t

    val gD = g.mapPairs {
      case ((i, s), 0.0) => mu_is_clipped(i, s)
      case ((i, s), 1.0) => 0.0
      case ((i, s), 2.0) => 1.0 - mu_is_clipped(i, s)
      case ((i, s), j) => throw new RuntimeException(s"Individual $i at snp $s has a genotype other than 0, 1, or 2: $j")
    }

    val k2 = for (i <- 0 until 4)
    yield for (j <- 0 until 4)
    yield {
      val numer = (for (s <- 0 until 8) yield
        (gD(i,s) - mu_is_clipped(i,s) * (1 - mu_is_clipped(i,s)) * (2.0 * phi_hat(i,i))) *
        (gD(j,s) - mu_is_clipped(j,s) * (1 - mu_is_clipped(j,s)) * (2.0 * phi_hat(j,j)))
      ).sum
      val denom = (for (s <- 0 until 8) yield
        mu_is_clipped(i, s) * (1.0 - mu_is_clipped(i, s)) *
        mu_is_clipped(j, s) * (1.0 - mu_is_clipped(j, s))
      ).sum
      numer / denom
    }

    println("k2")
    println(k2.map(_.mkString(", ")).mkString("\n"))
    println("")
  }
}


package secp256k1

import (
	"crypto/elliptic"
	"math/big"
	"sync"
)

type CurveParams struct {
	P       *big.Int
	N       *big.Int
	B       *big.Int
	Gx, Gy  *big.Int
	BitSize int
	Name    string
}

func (curve *CurveParams) Params() *elliptic.CurveParams {
	return &elliptic.CurveParams{
		P:       curve.P,
		N:       curve.N,
		B:       curve.B,
		Gx:      curve.Gx,
		Gy:      curve.Gy,
		BitSize: curve.BitSize,
		Name:    curve.Name,
	}
}

func (curve *CurveParams) IsOnCurve(x, y *big.Int) bool {
	// y^2
	y2 := new(big.Int).Mul(y, y)
	y2.Mod(y2, curve.P)

	// x^3 + b
	x3 := new(big.Int).Mul(x, x)
	x3.Mul(x3, x)
	x3.Add(x3, curve.B)
	x3.Mod(x3, curve.P)

	return y2.Cmp(x3) == 0
}

func zForAffine(xA, yA *big.Int) (z *big.Int) {
	z = new(big.Int)
	if xA.Sign() != 0 || yA.Sign() != 0 {
		return z.SetInt64(1)
	}
	return // point at infinity
}

func (curve *CurveParams) affineFromJacobian(x, y, z *big.Int) (xA, yA *big.Int) {
	// point at infinity
	if z.Sign() == 0 {
		return new(big.Int), new(big.Int)
	}

	// 1/Z^1
	zInv := new(big.Int).ModInverse(z, curve.P)

	// 1/Z^2
	zzInv := new(big.Int).Mul(zInv, zInv)

	// 1/Z^3
	zzzInv := new(big.Int).Mul(zzInv, zInv)

	// x = X/Z^2
	xA = new(big.Int).Mul(x, zzInv)
	xA.Mod(xA, curve.P)

	// y = Y/Z^3
	yA = new(big.Int).Mul(y, zzzInv)
	yA.Mod(yA, curve.P)

	return
}

func (curve *CurveParams) Add(x1, y1, x2, y2 *big.Int) (x, y *big.Int) {
	z1 := zForAffine(x1, y1)
	z2 := zForAffine(x2, y2)
	return curve.affineFromJacobian(curve.addJacobian(x1, y1, z1, x2, y2, z2))
}

// ref. https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-3.html#addition-add-2007-bl
func (curve *CurveParams) addJacobian(x1, y1, z1, x2, y2, z2 *big.Int) (x3, y3, z3 *big.Int) {
	x3, y3, z3 = new(big.Int), new(big.Int), new(big.Int)

	// (x1, y1) is point at infinity in affine coordinates
	if z1.Sign() == 0 {
		// (x3, y3) = (x2, y2) in affine coordinates
		x3.Set(x2)
		y3.Set(y2)
		z3.Set(z2)
		return
	}

	// (x2, y2) is point at infinity in affine coordinates
	if z2.Sign() == 0 {
		// (x3, y3) = (x1, y1) in affine coordinates
		x3.Set(x1)
		y3.Set(y1)
		z3.Set(z1)
		return
	}

	// Z1Z1 = Z1^2
	z1z1 := new(big.Int).Mul(z1, z1) // 1S

	// Z2Z2 = Z2^2
	z2z2 := new(big.Int).Mul(z2, z2) // 2S

	// U1 = X1*Z2Z2
	u1 := new(big.Int).Mul(x1, z2z2) // 1M
	u1.Mod(u1, curve.P)

	// U2 = X2*Z1Z1
	u2 := new(big.Int).Mul(x2, z1z1) // 2M
	u2.Mod(u2, curve.P)

	// S1 = Y1*Z2*Z2Z2
	s1 := new(big.Int).Mul(y1, z2) // 3M
	s1.Mul(s1, z2z2)               // 4M
	s1.Mod(s1, curve.P)

	// S2 = Y2*Z1*Z1Z1
	s2 := new(big.Int).Mul(y2, z1) // 5M
	s2.Mul(s2, z1z1)               // 6M
	s2.Mod(s2, curve.P)

	// x1 == x2 and y1 == y2 in affine coordinates
	if u1.Cmp(u2) == 0 && s1.Cmp(s2) == 0 {
		return curve.doubleJacobian(x1, y1, z1)
	}

	// H = U2 - U1
	h := new(big.Int).Sub(u2, u1) // 1add

	// I = (2*H)^2
	i := new(big.Int).Lsh(h, 1) // 1*2
	i.Mul(i, i)                 // 3S

	// J = H*I
	j := new(big.Int).Mul(h, i) // 7M

	// r = 2*(S2 - S1)
	r := new(big.Int).Sub(s2, s1) // 2add
	r.Lsh(r, 1)                   // 2*2

	// V = U1*I
	v := new(big.Int).Mul(u1, i) // 8M

	// tmp1 = 2*V
	tmp1 := new(big.Int).Lsh(v, 1) // 3*2

	// tmp2 = 2*S1*J
	tmp2 := new(big.Int).Mul(s1, j) // 9M
	tmp2.Lsh(tmp2, 1)               // 4*2

	// X3 = r^2 - J - 2*V
	x3.Mul(r, r)     // 4S
	x3.Sub(x3, j)    // 3add
	x3.Sub(x3, tmp1) // 4add
	x3.Mod(x3, curve.P)

	// Y3 = r*(V - X3) - 2*S1*J
	y3.Sub(v, x3)    // 5add
	y3.Mul(r, y3)    // 10M
	y3.Sub(y3, tmp2) // 6add
	y3.Mod(y3, curve.P)

	// Z3 = ((Z1 + Z2)^2 - Z1Z1 - Z2Z2)*H
	z3.Add(z1, z2)   // 7add
	z3.Mul(z3, z3)   // 5S
	z3.Sub(z3, z1z1) // 8add
	z3.Sub(z3, z2z2) // 9add
	z3.Mul(z3, h)    // 11M
	z3.Mod(z3, curve.P)

	// cost: 11M + 5S + 9add + 4*2

	return
}

func (curve *CurveParams) Double(x1, y1 *big.Int) (x, y *big.Int) {
	z1 := zForAffine(x1, y1)
	return curve.affineFromJacobian(curve.doubleJacobian(x1, y1, z1))
}

// ref. https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l
func (curve *CurveParams) doubleJacobian(x1, y1, z1 *big.Int) (x3, y3, z3 *big.Int) {
	x3, y3, z3 = new(big.Int), new(big.Int), new(big.Int)

	// (x1, y1) is point at infinity in affine coordinates
	if y1.Sign() == 0 || z1.Sign() == 0 {
		return
	}

	// A = X1^2
	a := new(big.Int).Mul(x1, x1) // 1S

	// B = Y1^2
	b := new(big.Int).Mul(y1, y1) // 2S

	// C = B^2
	c := new(big.Int).Mul(b, b) // 3S

	// D = 2*((X1 + B)^2 - A - C)
	d := new(big.Int).Add(x1, b) // 1add
	d.Mul(d, d)                  // 4S
	d.Sub(d, a)                  // 2add
	d.Sub(d, c)                  // 3add
	d.Lsh(d, 1)                  // 1*2

	// E = 3*A
	e := new(big.Int).Lsh(a, 1) // 2*2
	e.Add(e, a)                 // 4add

	// F = E^2
	f := new(big.Int).Mul(e, e) // 5S

	// tmp1 = 2*D
	tmp1 := new(big.Int).Lsh(d, 1) // 3*2

	// tmp2 = 8*C
	tmp2 := new(big.Int).Lsh(c, 3) // 1*8

	// X3 = F - 2*D
	x3.Sub(f, tmp1) // 5add
	x3.Mod(x3, curve.P)

	// Y3 = E * (D - X3) - 8*C
	y3.Sub(d, x3)    // 6add
	y3.Mul(e, y3)    // 1M
	y3.Sub(y3, tmp2) // 7add
	y3.Mod(y3, curve.P)

	// Z3 = 2 * Y1 * Z1
	z3.Mul(y1, z1) // 2M
	z3.Lsh(z3, 1)  // 4*2
	z3.Mod(z3, curve.P)

	// cost: 2M + 5S + 7add + 4*2 + 1*8

	return
}

func (curve *CurveParams) ScalarMult(x1, y1 *big.Int, k []byte) (x, y *big.Int) {
	z1 := new(big.Int).SetInt64(1)
	x, y, z := new(big.Int), new(big.Int), new(big.Int)

	for _, byte := range k {
		for bitNum := 0; bitNum < 8; bitNum++ {
			x, y, z = curve.doubleJacobian(x, y, z)
			if byte&0x80 == 0x80 {
				x, y, z = curve.addJacobian(x1, y1, z1, x, y, z)
			}
			byte <<= 1
		}
	}

	return curve.affineFromJacobian(x, y, z)
}

func (curve *CurveParams) ScalarBaseMult(k []byte) (x, y *big.Int) {
	return curve.ScalarMult(curve.Gx, curve.Gy, k)
}

const (
	p       = "fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f"
	n       = "fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141"
	b       = "0000000000000000000000000000000000000000000000000000000000000007"
	gx      = "79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798"
	gy      = "483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8"
	bitSize = 256
)

var (
	initonce  sync.Once
	secp256k1 *CurveParams
)

func initS256() {
	secp256k1 = &CurveParams{Name: "secp256k1"}
	secp256k1.P = fromHexToBigInt(p)
	secp256k1.N = fromHexToBigInt(n)
	secp256k1.B = fromHexToBigInt(b)
	secp256k1.Gx = fromHexToBigInt(gx)
	secp256k1.Gy = fromHexToBigInt(gy)
	secp256k1.BitSize = bitSize
}

func S256() elliptic.Curve {
	initonce.Do(initS256)
	return secp256k1
}

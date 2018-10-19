package secp256k1

import "math/big"

func hexToBigInt(s string) *big.Int {
	i, ok := new(big.Int).SetString(s, 16)
	if !ok {
		panic("invalid hex : " + s)
	}
	return i
}

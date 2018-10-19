package secp256k1

import (
	"math/big"
	"testing"
)

func TestHexToBigInt(t *testing.T) {
	expected := new(big.Int).SetInt64(81985529216486895)

	actual := hexToBigInt("0123456789abcdef")
	if actual.Cmp(expected) != 0 {
		t.Errorf("expected: %x, actual: %x", expected, actual)
	}
}

func TestHexToBigInt_panic(t *testing.T) {
	defer func() {
		if r := recover(); r == nil {
			t.Errorf("no panic")
		}
	}()

	hexToBigInt("0123456789abcdefg")
}

package secp256k1

import (
	"math/big"
	"testing"
)

func TestFromHexToBigInt(t *testing.T) {
	expected := new(big.Int).SetInt64(81985529216486895)

	actual := fromHexToBigInt("0123456789abcdef")
	if actual.Cmp(expected) != 0 {
		t.Errorf("expected: %x, actual: %x", expected, actual)
	}
}

func TestFromHexToBigInt_panic(t *testing.T) {
	defer func() {
		if r := recover(); r == nil {
			t.Errorf("no panic")
		}
	}()

	fromHexToBigInt("0123456789abcdefg")
}

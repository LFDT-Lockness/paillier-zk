# paillier-zk

This crate provides ZK-proofs for some properties about paillier encryption.
See the module docs for the properties and examples of usage.

This library is built on top of [fast-paillier](https://lib.rs/fast-paillier) crate.
This crate and the underlying big integer implementation are reexported for the
consumer to be able to use them, instead of trying to match a version.

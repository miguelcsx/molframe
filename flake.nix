{
  description = "MolFrame development environment";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-26.05";

    rust-overlay = {
      url = "github:oxalica/rust-overlay";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs =
    {
      nixpkgs,
      rust-overlay,
      ...
    }:
    let
      systems = [
        "aarch64-darwin"
        "x86_64-darwin"
        "aarch64-linux"
        "x86_64-linux"
      ];

      forAllSystems = nixpkgs.lib.genAttrs systems;
    in
    {
      devShells = forAllSystems (
        system:
        let
          pkgs = import nixpkgs {
            inherit system;
            overlays = [ rust-overlay.overlays.default ];
          };

          rust = pkgs.rust-bin.fromRustupToolchainFile ./rust-toolchain.toml;
          python = pkgs.python313;
        in
        {
          default = pkgs.mkShell {
            packages = [
              rust

              # Python / PyO3 toolchain
              python
              pkgs.uv
              pkgs.maturin

              # Repository / native build utilities
              pkgs.git
              pkgs.pkg-config
              pkgs.jq
            ];

            # Keep uv on the same interpreter as the repository contract.
            UV_PYTHON = "${python}/bin/python3";

            RUST_BACKTRACE = "1";
            CARGO_TERM_COLOR = "always";
          };
        }
      );
    };
}

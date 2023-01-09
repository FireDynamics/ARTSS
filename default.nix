{pkgs ? import <nixpkgs> {}}:

let 
  fds = p: with p; [
    (
      buildPythonPackage rec {
        pname = "fdsreader";
        version = "1.9.3";
        src = pkgs.fetchFromGitHub {
          owner = "FireDynamics";
          repo = "${name}";
          rev = "${version}";
          sha256 = "1b3434pd95ynd441a6asd1666808293m0dfdgxgv4sc8hybsq973";
        };

        propagatedBuildInputs = with p; [
          numpy
        ];
        doCheck = false;
      }
    )
  ];

  myPythonPackages = ps: with ps; [
    flake8
    pylint
    pydantic
    python-lsp-server
    numpy
    scipy
    pandas
    h5py
    retry
    incremental
    matplotlib
    tqdm
    fds
  ];

in

  pkgs.mkShell {
    buildInputs = [
      (pkgs.python3.withPackages myPythonPackages)
      pkgs.ccls
      pkgs.cmake
      pkgs.ctags
      pkgs.gcc
      pkgs.gdb
      pkgs.git
      pkgs.gitAndTools.gh
      pkgs.hdf5
      pkgs.jq
      pkgs.llvmPackages.bintools
      pkgs.openmpi
      pkgs.libzip
      pkgs.nvim-dev
      pkgs.pkg-config 
      pkgs.pkgconfig
    ];
}

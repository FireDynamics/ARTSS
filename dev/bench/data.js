window.BENCHMARK_DATA = {
  "lastUpdate": 1652020450675,
  "repoUrl": "https://github.com/FireDynamics/ARTSS",
  "entries": {
    "Benchmark": [
      {
        "commit": {
          "author": {
            "email": "c.von.mach@fz-juelich.de",
            "name": "Christian von Mach",
            "username": "ChristianVonMach"
          },
          "committer": {
            "email": "c.von.mach@fz-juelich.de",
            "name": "Christian von Mach",
            "username": "ChristianVonMach"
          },
          "distinct": true,
          "id": "283bfd4d70b328c840e69271f93de156d5c563eb",
          "message": "fixed double output-file",
          "timestamp": "2022-05-08T16:21:59+02:00",
          "tree_id": "243db0f8229295fe0ff27a6169fc6628b42fe94f",
          "url": "https://github.com/FireDynamics/ARTSS/commit/283bfd4d70b328c840e69271f93de156d5c563eb"
        },
        "date": 1652020450081,
        "tool": "googlecpp",
        "benches": [
          {
            "name": "BM_AddScalar/8",
            "value": 1128.1626946323634,
            "unit": "ns/iter",
            "extra": "iterations: 637163\ncpu: 1127.6343729940377 ns\nthreads: 1"
          },
          {
            "name": "BM_AddScalar/64",
            "value": 1156.4804227151226,
            "unit": "ns/iter",
            "extra": "iterations: 606366\ncpu: 1156.4329134549098 ns\nthreads: 1"
          },
          {
            "name": "BM_AddScalar/512",
            "value": 1381.860972628759,
            "unit": "ns/iter",
            "extra": "iterations: 506483\ncpu: 1381.707974403879 ns\nthreads: 1"
          },
          {
            "name": "BM_AddScalar/4096",
            "value": 3127.3440244436,
            "unit": "ns/iter",
            "extra": "iterations: 223862\ncpu: 3127.1354673861574 ns\nthreads: 1"
          },
          {
            "name": "BM_AddScalar/32768",
            "value": 15779.68152651082,
            "unit": "ns/iter",
            "extra": "iterations: 44415\ncpu: 15778.655859506913 ns\nthreads: 1"
          },
          {
            "name": "BM_AddScalar/262144",
            "value": 168920.511470675,
            "unit": "ns/iter",
            "extra": "iterations: 4141\ncpu: 168901.4730741367 ns\nthreads: 1"
          },
          {
            "name": "BM_AddScalar/2097152",
            "value": 1361383.9319066613,
            "unit": "ns/iter",
            "extra": "iterations: 514\ncpu: 1361171.4007782086 ns\nthreads: 1"
          },
          {
            "name": "BM_AddScalar/8388608",
            "value": 23632611.83333331,
            "unit": "ns/iter",
            "extra": "iterations: 30\ncpu: 23628436.66666666 ns\nthreads: 1"
          },
          {
            "name": "BM_AddFields/8",
            "value": 2280.2078182700448,
            "unit": "ns/iter",
            "extra": "iterations: 305464\ncpu: 2279.17528743158 ns\nthreads: 1"
          },
          {
            "name": "BM_AddFields/64",
            "value": 2328.127126174572,
            "unit": "ns/iter",
            "extra": "iterations: 302715\ncpu: 2316.854467072992 ns\nthreads: 1"
          },
          {
            "name": "BM_AddFields/512",
            "value": 2679.2414147004297,
            "unit": "ns/iter",
            "extra": "iterations: 262105\ncpu: 2677.31214589573 ns\nthreads: 1"
          },
          {
            "name": "BM_AddFields/4096",
            "value": 5718.390839072327,
            "unit": "ns/iter",
            "extra": "iterations: 122695\ncpu: 5716.902074249155 ns\nthreads: 1"
          },
          {
            "name": "BM_AddFields/32768",
            "value": 28399.272118505236,
            "unit": "ns/iter",
            "extra": "iterations: 24640\ncpu: 28396.724837662343 ns\nthreads: 1"
          },
          {
            "name": "BM_AddFields/262144",
            "value": 445741.94652404037,
            "unit": "ns/iter",
            "extra": "iterations: 2057\ncpu: 347474.4773942636 ns\nthreads: 1"
          },
          {
            "name": "BM_AddFields/2097152",
            "value": 9454695.298701722,
            "unit": "ns/iter",
            "extra": "iterations: 77\ncpu: 9453155.844155842 ns\nthreads: 1"
          },
          {
            "name": "BM_AddFields/8388608",
            "value": 84322609.37499337,
            "unit": "ns/iter",
            "extra": "iterations: 8\ncpu: 84315875.00000015 ns\nthreads: 1"
          },
          {
            "name": "BM_MulScalar/8",
            "value": 1136.3937955456547,
            "unit": "ns/iter",
            "extra": "iterations: 614494\ncpu: 1136.3036905160986 ns\nthreads: 1"
          },
          {
            "name": "BM_MulScalar/64",
            "value": 1165.0539367969875,
            "unit": "ns/iter",
            "extra": "iterations: 600666\ncpu: 1164.9310931532664 ns\nthreads: 1"
          },
          {
            "name": "BM_MulScalar/512",
            "value": 1388.2855578180574,
            "unit": "ns/iter",
            "extra": "iterations: 505810\ncpu: 1388.1955675055872 ns\nthreads: 1"
          },
          {
            "name": "BM_MulScalar/4096",
            "value": 3115.201465201406,
            "unit": "ns/iter",
            "extra": "iterations: 224952\ncpu: 3114.930296240969 ns\nthreads: 1"
          },
          {
            "name": "BM_MulScalar/32768",
            "value": 15782.174564435383,
            "unit": "ns/iter",
            "extra": "iterations: 44253\ncpu: 15781.067950195464 ns\nthreads: 1"
          },
          {
            "name": "BM_MulScalar/262144",
            "value": 179556.80353938352,
            "unit": "ns/iter",
            "extra": "iterations: 3899\ncpu: 179543.8830469348 ns\nthreads: 1"
          },
          {
            "name": "BM_MulScalar/2097152",
            "value": 1422902.2525458096,
            "unit": "ns/iter",
            "extra": "iterations: 491\ncpu: 1422872.5050916476 ns\nthreads: 1"
          },
          {
            "name": "BM_MulScalar/8388608",
            "value": 45535981.26666808,
            "unit": "ns/iter",
            "extra": "iterations: 15\ncpu: 45530859.999999985 ns\nthreads: 1"
          },
          {
            "name": "BM_MulFields/8",
            "value": 2269.5657913096015,
            "unit": "ns/iter",
            "extra": "iterations: 308179\ncpu: 2269.0549972580798 ns\nthreads: 1"
          },
          {
            "name": "BM_MulFields/64",
            "value": 2350.7101570252157,
            "unit": "ns/iter",
            "extra": "iterations: 298296\ncpu: 2350.0908493576944 ns\nthreads: 1"
          },
          {
            "name": "BM_MulFields/512",
            "value": 3022.44534369222,
            "unit": "ns/iter",
            "extra": "iterations: 231675\ncpu: 3021.483543757404 ns\nthreads: 1"
          },
          {
            "name": "BM_MulFields/4096",
            "value": 8522.305693550745,
            "unit": "ns/iter",
            "extra": "iterations: 82128\ncpu: 8521.241233196954 ns\nthreads: 1"
          },
          {
            "name": "BM_MulFields/32768",
            "value": 50643.386886672,
            "unit": "ns/iter",
            "extra": "iterations: 13818\ncpu: 50630.95238095235 ns\nthreads: 1"
          },
          {
            "name": "BM_MulFields/262144",
            "value": 534358.105945192,
            "unit": "ns/iter",
            "extra": "iterations: 1312\ncpu: 534289.024390244 ns\nthreads: 1"
          },
          {
            "name": "BM_MulFields/2097152",
            "value": 16578243.67499927,
            "unit": "ns/iter",
            "extra": "iterations: 40\ncpu: 16576422.499999931 ns\nthreads: 1"
          },
          {
            "name": "BM_MulFields/8388608",
            "value": 103492455.42857131,
            "unit": "ns/iter",
            "extra": "iterations: 7\ncpu: 103485214.28571442 ns\nthreads: 1"
          }
        ]
      }
    ]
  }
}
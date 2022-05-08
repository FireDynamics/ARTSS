window.BENCHMARK_DATA = {
  "lastUpdate": 1652021137178,
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
      },
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
          "id": "05b0b6fcefa53027e472af8811f4ad5456a7ddf1",
          "message": "removed dummy loop to Field mul",
          "timestamp": "2022-05-08T16:37:47+02:00",
          "tree_id": "e32e736f971768c4f314c760cf28f1c9133ba29e",
          "url": "https://github.com/FireDynamics/ARTSS/commit/05b0b6fcefa53027e472af8811f4ad5456a7ddf1"
        },
        "date": 1652021136449,
        "tool": "googlecpp",
        "benches": [
          {
            "name": "BM_AddScalar/8",
            "value": 1132.9074361775242,
            "unit": "ns/iter",
            "extra": "iterations: 641741\ncpu: 1132.7660224296096 ns\nthreads: 1"
          },
          {
            "name": "BM_AddScalar/64",
            "value": 1161.3212907350792,
            "unit": "ns/iter",
            "extra": "iterations: 603718\ncpu: 1161.257408260148 ns\nthreads: 1"
          },
          {
            "name": "BM_AddScalar/512",
            "value": 1384.2283134448921,
            "unit": "ns/iter",
            "extra": "iterations: 505901\ncpu: 1384.0925398447523 ns\nthreads: 1"
          },
          {
            "name": "BM_AddScalar/4096",
            "value": 3139.6388952465168,
            "unit": "ns/iter",
            "extra": "iterations: 222819\ncpu: 3138.886719714208 ns\nthreads: 1"
          },
          {
            "name": "BM_AddScalar/32768",
            "value": 15523.551897272528,
            "unit": "ns/iter",
            "extra": "iterations: 45012\ncpu: 15513.072069670308 ns\nthreads: 1"
          },
          {
            "name": "BM_AddScalar/262144",
            "value": 174458.92924763376,
            "unit": "ns/iter",
            "extra": "iterations: 4014\ncpu: 174452.9147982061 ns\nthreads: 1"
          },
          {
            "name": "BM_AddScalar/2097152",
            "value": 1495664.554140111,
            "unit": "ns/iter",
            "extra": "iterations: 471\ncpu: 1495550.3184713374 ns\nthreads: 1"
          },
          {
            "name": "BM_AddScalar/8388608",
            "value": 26869796.480000332,
            "unit": "ns/iter",
            "extra": "iterations: 25\ncpu: 26868200.00000001 ns\nthreads: 1"
          },
          {
            "name": "BM_AddFields/8",
            "value": 2278.728652224797,
            "unit": "ns/iter",
            "extra": "iterations: 306238\ncpu: 2277.801252620512 ns\nthreads: 1"
          },
          {
            "name": "BM_AddFields/64",
            "value": 2319.1898506572447,
            "unit": "ns/iter",
            "extra": "iterations: 301990\ncpu: 2318.0072187820806 ns\nthreads: 1"
          },
          {
            "name": "BM_AddFields/512",
            "value": 2682.043499003712,
            "unit": "ns/iter",
            "extra": "iterations: 262006\ncpu: 2681.2336358709326 ns\nthreads: 1"
          },
          {
            "name": "BM_AddFields/4096",
            "value": 5664.646273385107,
            "unit": "ns/iter",
            "extra": "iterations: 123423\ncpu: 5663.458998727954 ns\nthreads: 1"
          },
          {
            "name": "BM_AddFields/32768",
            "value": 28155.073410383586,
            "unit": "ns/iter",
            "extra": "iterations: 24833\ncpu: 28152.33761526999 ns\nthreads: 1"
          },
          {
            "name": "BM_AddFields/262144",
            "value": 344206.5667976157,
            "unit": "ns/iter",
            "extra": "iterations: 2036\ncpu: 344168.61493123707 ns\nthreads: 1"
          },
          {
            "name": "BM_AddFields/2097152",
            "value": 10565033.712122262,
            "unit": "ns/iter",
            "extra": "iterations: 66\ncpu: 10564319.696969697 ns\nthreads: 1"
          },
          {
            "name": "BM_AddFields/8388608",
            "value": 79496276.33334129,
            "unit": "ns/iter",
            "extra": "iterations: 9\ncpu: 79484433.33333345 ns\nthreads: 1"
          },
          {
            "name": "BM_MulScalar/8",
            "value": 1137.949896429136,
            "unit": "ns/iter",
            "extra": "iterations: 613589\ncpu: 1137.9165858579606 ns\nthreads: 1"
          },
          {
            "name": "BM_MulScalar/64",
            "value": 1166.40366941938,
            "unit": "ns/iter",
            "extra": "iterations: 600858\ncpu: 1166.2151456750212 ns\nthreads: 1"
          },
          {
            "name": "BM_MulScalar/512",
            "value": 1382.030415871045,
            "unit": "ns/iter",
            "extra": "iterations: 504868\ncpu: 1381.9390414920322 ns\nthreads: 1"
          },
          {
            "name": "BM_MulScalar/4096",
            "value": 3144.2487056885816,
            "unit": "ns/iter",
            "extra": "iterations: 222319\ncpu: 3143.99803885407 ns\nthreads: 1"
          },
          {
            "name": "BM_MulScalar/32768",
            "value": 15469.10367930112,
            "unit": "ns/iter",
            "extra": "iterations: 45226\ncpu: 15468.274886127407 ns\nthreads: 1"
          },
          {
            "name": "BM_MulScalar/262144",
            "value": 175158.0797912027,
            "unit": "ns/iter",
            "extra": "iterations: 4023\ncpu: 174508.55083271195 ns\nthreads: 1"
          },
          {
            "name": "BM_MulScalar/2097152",
            "value": 1394322.0732674282,
            "unit": "ns/iter",
            "extra": "iterations: 505\ncpu: 1394090.8910891123 ns\nthreads: 1"
          },
          {
            "name": "BM_MulScalar/8388608",
            "value": 43287255.18750076,
            "unit": "ns/iter",
            "extra": "iterations: 16\ncpu: 43283793.750000134 ns\nthreads: 1"
          },
          {
            "name": "BM_MulFields/8",
            "value": 2282.22150216776,
            "unit": "ns/iter",
            "extra": "iterations: 306304\ncpu: 2281.658091307982 ns\nthreads: 1"
          },
          {
            "name": "BM_MulFields/64",
            "value": 2308.595496827247,
            "unit": "ns/iter",
            "extra": "iterations: 300366\ncpu: 2307.9403128183662 ns\nthreads: 1"
          },
          {
            "name": "BM_MulFields/512",
            "value": 2675.7648178303816,
            "unit": "ns/iter",
            "extra": "iterations: 262201\ncpu: 2674.8387687308627 ns\nthreads: 1"
          },
          {
            "name": "BM_MulFields/4096",
            "value": 5665.218839223222,
            "unit": "ns/iter",
            "extra": "iterations: 123986\ncpu: 5663.562821608879 ns\nthreads: 1"
          },
          {
            "name": "BM_MulFields/32768",
            "value": 28109.38629596808,
            "unit": "ns/iter",
            "extra": "iterations: 24854\ncpu: 28104.864408143552 ns\nthreads: 1"
          },
          {
            "name": "BM_MulFields/262144",
            "value": 344511.2693254513,
            "unit": "ns/iter",
            "extra": "iterations: 2031\ncpu: 344464.10635155055 ns\nthreads: 1"
          },
          {
            "name": "BM_MulFields/2097152",
            "value": 11544502.075471943,
            "unit": "ns/iter",
            "extra": "iterations: 53\ncpu: 11539332.075471707 ns\nthreads: 1"
          },
          {
            "name": "BM_MulFields/8388608",
            "value": 82289826.33334504,
            "unit": "ns/iter",
            "extra": "iterations: 9\ncpu: 82252677.77777785 ns\nthreads: 1"
          }
        ]
      }
    ]
  }
}
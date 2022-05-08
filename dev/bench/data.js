window.BENCHMARK_DATA = {
  "lastUpdate": 1652022825200,
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
          "id": "abc78455c44076785f50863f00ddd0970c60c4a4",
          "message": "multiple threads",
          "timestamp": "2022-05-08T17:02:18+02:00",
          "tree_id": "29cd87596b0ad04e378058968a859e35fb996ea9",
          "url": "https://github.com/FireDynamics/ARTSS/commit/abc78455c44076785f50863f00ddd0970c60c4a4"
        },
        "date": 1652022823240,
        "tool": "googlecpp",
        "benches": [
          {
            "name": "BM_AddScalar/8/threads:1",
            "value": 1533.3670667083302,
            "unit": "ns/iter",
            "extra": "iterations: 476867\ncpu: 1533.239666405937 ns\nthreads: 1"
          },
          {
            "name": "BM_AddScalar/8/threads:2",
            "value": 2210.023815043873,
            "unit": "ns/iter",
            "extra": "iterations: 194100\ncpu: 3869.460587326121 ns\nthreads: 2"
          },
          {
            "name": "BM_AddScalar/8/threads:4",
            "value": 2123.9978743505394,
            "unit": "ns/iter",
            "extra": "iterations: 186296\ncpu: 3717.032571821189 ns\nthreads: 4"
          },
          {
            "name": "BM_AddScalar/8/threads:8",
            "value": 2004.2139061406926,
            "unit": "ns/iter",
            "extra": "iterations: 182912\ncpu: 3935.6280615815263 ns\nthreads: 8"
          },
          {
            "name": "BM_AddScalar/64/threads:1",
            "value": 1505.3049825573596,
            "unit": "ns/iter",
            "extra": "iterations: 462092\ncpu: 1505.1249534724684 ns\nthreads: 1"
          },
          {
            "name": "BM_AddScalar/64/threads:2",
            "value": 2218.3963635592345,
            "unit": "ns/iter",
            "extra": "iterations: 188206\ncpu: 3863.0660021465833 ns\nthreads: 2"
          },
          {
            "name": "BM_AddScalar/64/threads:4",
            "value": 2055.2621088063215,
            "unit": "ns/iter",
            "extra": "iterations: 329760\ncpu: 3663.352438136827 ns\nthreads: 4"
          },
          {
            "name": "BM_AddScalar/64/threads:8",
            "value": 2178.7799123282393,
            "unit": "ns/iter",
            "extra": "iterations: 165960\ncpu: 4126.757652446374 ns\nthreads: 8"
          },
          {
            "name": "BM_AddScalar/512/threads:1",
            "value": 1922.9948358460947,
            "unit": "ns/iter",
            "extra": "iterations: 383025\ncpu: 1921.1884341753155 ns\nthreads: 1"
          },
          {
            "name": "BM_AddScalar/512/threads:2",
            "value": 2367.8932711709167,
            "unit": "ns/iter",
            "extra": "iterations: 167518\ncpu: 4174.330519705344 ns\nthreads: 2"
          },
          {
            "name": "BM_AddScalar/512/threads:4",
            "value": 2170.0237920841078,
            "unit": "ns/iter",
            "extra": "iterations: 184988\ncpu: 3532.244794256921 ns\nthreads: 4"
          },
          {
            "name": "BM_AddScalar/512/threads:8",
            "value": 2095.0937929041015,
            "unit": "ns/iter",
            "extra": "iterations: 174808\ncpu: 4025.914145805684 ns\nthreads: 8"
          },
          {
            "name": "BM_AddScalar/4096/threads:1",
            "value": 4173.992845702655,
            "unit": "ns/iter",
            "extra": "iterations: 169828\ncpu: 4173.605059236405 ns\nthreads: 1"
          },
          {
            "name": "BM_AddScalar/4096/threads:2",
            "value": 3451.8932671653542,
            "unit": "ns/iter",
            "extra": "iterations: 108008\ncpu: 6302.915524775945 ns\nthreads: 2"
          },
          {
            "name": "BM_AddScalar/4096/threads:4",
            "value": 3510.9793345312273,
            "unit": "ns/iter",
            "extra": "iterations: 110208\ncpu: 6474.237804878048 ns\nthreads: 4"
          },
          {
            "name": "BM_AddScalar/4096/threads:8",
            "value": 3368.794876426295,
            "unit": "ns/iter",
            "extra": "iterations: 107264\ncpu: 6095.2537664081165 ns\nthreads: 8"
          },
          {
            "name": "BM_AddScalar/32768/threads:1",
            "value": 27832.250960705827,
            "unit": "ns/iter",
            "extra": "iterations: 25502\ncpu: 27787.816641831996 ns\nthreads: 1"
          },
          {
            "name": "BM_AddScalar/32768/threads:2",
            "value": 15081.376186202877,
            "unit": "ns/iter",
            "extra": "iterations: 23394\ncpu: 28998.08497905449 ns\nthreads: 2"
          },
          {
            "name": "BM_AddScalar/32768/threads:4",
            "value": 14644.1173704169,
            "unit": "ns/iter",
            "extra": "iterations: 24772\ncpu: 28343.64605199422 ns\nthreads: 4"
          },
          {
            "name": "BM_AddScalar/32768/threads:8",
            "value": 14039.37725752506,
            "unit": "ns/iter",
            "extra": "iterations: 23920\ncpu: 28932.679765886292 ns\nthreads: 8"
          },
          {
            "name": "BM_AddScalar/262144/threads:1",
            "value": 263072.8185255111,
            "unit": "ns/iter",
            "extra": "iterations: 2645\ncpu: 263051.6446124766 ns\nthreads: 1"
          },
          {
            "name": "BM_AddScalar/262144/threads:2",
            "value": 132287.63135279043,
            "unit": "ns/iter",
            "extra": "iterations: 2676\ncpu: 264142.0777279519 ns\nthreads: 2"
          },
          {
            "name": "BM_AddScalar/262144/threads:4",
            "value": 137355.03687889277,
            "unit": "ns/iter",
            "extra": "iterations: 2576\ncpu: 283400.85403726716 ns\nthreads: 4"
          },
          {
            "name": "BM_AddScalar/262144/threads:8",
            "value": 131515.23621616478,
            "unit": "ns/iter",
            "extra": "iterations: 2512\ncpu: 276625.9952229302 ns\nthreads: 8"
          },
          {
            "name": "BM_AddScalar/2097152/threads:1",
            "value": 2355911.310000162,
            "unit": "ns/iter",
            "extra": "iterations: 300\ncpu: 2355271.333333331 ns\nthreads: 1"
          },
          {
            "name": "BM_AddScalar/2097152/threads:2",
            "value": 1322425.2520001302,
            "unit": "ns/iter",
            "extra": "iterations: 250\ncpu: 2640093.5999999987 ns\nthreads: 2"
          },
          {
            "name": "BM_AddScalar/2097152/threads:4",
            "value": 1738691.17219399,
            "unit": "ns/iter",
            "extra": "iterations: 196\ncpu: 3567812.7551020426 ns\nthreads: 4"
          },
          {
            "name": "BM_AddScalar/2097152/threads:8",
            "value": 1745490.2492559655,
            "unit": "ns/iter",
            "extra": "iterations: 168\ncpu: 4062098.809523807 ns\nthreads: 8"
          },
          {
            "name": "BM_AddScalar/8388608/threads:1",
            "value": 57087521.384617575,
            "unit": "ns/iter",
            "extra": "iterations: 13\ncpu: 57031992.30769229 ns\nthreads: 1"
          },
          {
            "name": "BM_AddScalar/8388608/threads:2",
            "value": 30211139.16666707,
            "unit": "ns/iter",
            "extra": "iterations: 12\ncpu: 60037416.66666674 ns\nthreads: 2"
          },
          {
            "name": "BM_AddScalar/8388608/threads:4",
            "value": 28466618.854163054,
            "unit": "ns/iter",
            "extra": "iterations: 12\ncpu: 58341558.33333325 ns\nthreads: 4"
          },
          {
            "name": "BM_AddScalar/8388608/threads:8",
            "value": 28981824.05468841,
            "unit": "ns/iter",
            "extra": "iterations: 16\ncpu: 58915231.249999985 ns\nthreads: 8"
          },
          {
            "name": "BM_AddFields/8/threads:1",
            "value": 2935.5206640337465,
            "unit": "ns/iter",
            "extra": "iterations: 238482\ncpu: 2935.376254811691 ns\nthreads: 1"
          },
          {
            "name": "BM_AddFields/8/threads:2",
            "value": 4438.682179612856,
            "unit": "ns/iter",
            "extra": "iterations: 85896\ncpu: 7768.992735400952 ns\nthreads: 2"
          },
          {
            "name": "BM_AddFields/8/threads:4",
            "value": 4036.8942954776485,
            "unit": "ns/iter",
            "extra": "iterations: 104908\ncpu: 7297.323369047157 ns\nthreads: 4"
          },
          {
            "name": "BM_AddFields/8/threads:8",
            "value": 4202.282223024716,
            "unit": "ns/iter",
            "extra": "iterations: 81704\ncpu: 8618.404239694511 ns\nthreads: 8"
          },
          {
            "name": "BM_AddFields/64/threads:1",
            "value": 2946.8390974303297,
            "unit": "ns/iter",
            "extra": "iterations: 239760\ncpu: 2946.415999332661 ns\nthreads: 1"
          },
          {
            "name": "BM_AddFields/64/threads:2",
            "value": 4508.392533363772,
            "unit": "ns/iter",
            "extra": "iterations: 86920\ncpu: 7780.283018867912 ns\nthreads: 2"
          },
          {
            "name": "BM_AddFields/64/threads:4",
            "value": 4181.637024141827,
            "unit": "ns/iter",
            "extra": "iterations: 95932\ncpu: 7158.420547888109 ns\nthreads: 4"
          },
          {
            "name": "BM_AddFields/64/threads:8",
            "value": 4227.64726185498,
            "unit": "ns/iter",
            "extra": "iterations: 81232\ncpu: 8322.19937955483 ns\nthreads: 8"
          },
          {
            "name": "BM_AddFields/512/threads:1",
            "value": 3694.713064674299,
            "unit": "ns/iter",
            "extra": "iterations: 190353\ncpu: 3693.700650895963 ns\nthreads: 1"
          },
          {
            "name": "BM_AddFields/512/threads:2",
            "value": 4840.5912841141035,
            "unit": "ns/iter",
            "extra": "iterations: 82298\ncpu: 8449.019417239795 ns\nthreads: 2"
          },
          {
            "name": "BM_AddFields/512/threads:4",
            "value": 4687.765112135637,
            "unit": "ns/iter",
            "extra": "iterations: 95420\ncpu: 8214.961224062052 ns\nthreads: 4"
          },
          {
            "name": "BM_AddFields/512/threads:8",
            "value": 4526.329604496279,
            "unit": "ns/iter",
            "extra": "iterations: 82376\ncpu: 8462.706370787613 ns\nthreads: 8"
          },
          {
            "name": "BM_AddFields/4096/threads:1",
            "value": 8509.888800972672,
            "unit": "ns/iter",
            "extra": "iterations: 84677\ncpu: 8508.056497041705 ns\nthreads: 1"
          },
          {
            "name": "BM_AddFields/4096/threads:2",
            "value": 6989.390517178602,
            "unit": "ns/iter",
            "extra": "iterations: 55068\ncpu: 12513.492409384791 ns\nthreads: 2"
          },
          {
            "name": "BM_AddFields/4096/threads:4",
            "value": 7228.480687500393,
            "unit": "ns/iter",
            "extra": "iterations: 40000\ncpu: 13117.132499999996 ns\nthreads: 4"
          },
          {
            "name": "BM_AddFields/4096/threads:8",
            "value": 6568.049585816644,
            "unit": "ns/iter",
            "extra": "iterations: 57040\ncpu: 12467.945301542803 ns\nthreads: 8"
          },
          {
            "name": "BM_AddFields/32768/threads:1",
            "value": 66862.80024969141,
            "unit": "ns/iter",
            "extra": "iterations: 10413\ncpu: 66789.64755593966 ns\nthreads: 1"
          },
          {
            "name": "BM_AddFields/32768/threads:2",
            "value": 32954.093011417506,
            "unit": "ns/iter",
            "extra": "iterations: 9816\ncpu: 64900.79462102686 ns\nthreads: 2"
          },
          {
            "name": "BM_AddFields/32768/threads:4",
            "value": 31297.279535677153,
            "unit": "ns/iter",
            "extra": "iterations: 11156\ncpu: 63420.50914306188 ns\nthreads: 4"
          },
          {
            "name": "BM_AddFields/32768/threads:8",
            "value": 31922.43367187331,
            "unit": "ns/iter",
            "extra": "iterations: 8000\ncpu: 65693.17499999986 ns\nthreads: 8"
          },
          {
            "name": "BM_AddFields/262144/threads:1",
            "value": 484632.29341742664,
            "unit": "ns/iter",
            "extra": "iterations: 1428\ncpu: 484545.3781512591 ns\nthreads: 1"
          },
          {
            "name": "BM_AddFields/262144/threads:2",
            "value": 241576.00916231825,
            "unit": "ns/iter",
            "extra": "iterations: 1528\ncpu: 481890.7068062826 ns\nthreads: 2"
          },
          {
            "name": "BM_AddFields/262144/threads:4",
            "value": 234775.74914851054,
            "unit": "ns/iter",
            "extra": "iterations: 1468\ncpu: 478530.8583106267 ns\nthreads: 4"
          },
          {
            "name": "BM_AddFields/262144/threads:8",
            "value": 247478.2407507133,
            "unit": "ns/iter",
            "extra": "iterations: 1392\ncpu: 512784.7701149438 ns\nthreads: 8"
          },
          {
            "name": "BM_AddFields/2097152/threads:1",
            "value": 16170664.933333177,
            "unit": "ns/iter",
            "extra": "iterations: 45\ncpu: 16156575.555555562 ns\nthreads: 1"
          },
          {
            "name": "BM_AddFields/2097152/threads:2",
            "value": 6555504.107142595,
            "unit": "ns/iter",
            "extra": "iterations: 56\ncpu: 13074558.928571463 ns\nthreads: 2"
          },
          {
            "name": "BM_AddFields/2097152/threads:4",
            "value": 5190397.379464632,
            "unit": "ns/iter",
            "extra": "iterations: 56\ncpu: 12126130.357142849 ns\nthreads: 4"
          },
          {
            "name": "BM_AddFields/2097152/threads:8",
            "value": 4958416.730468374,
            "unit": "ns/iter",
            "extra": "iterations: 64\ncpu: 11831806.25000004 ns\nthreads: 8"
          },
          {
            "name": "BM_AddFields/8388608/threads:1",
            "value": 102214808.49999882,
            "unit": "ns/iter",
            "extra": "iterations: 6\ncpu: 102200483.33333291 ns\nthreads: 1"
          },
          {
            "name": "BM_AddFields/8388608/threads:2",
            "value": 54241531.41666466,
            "unit": "ns/iter",
            "extra": "iterations: 6\ncpu: 107020950.00000042 ns\nthreads: 2"
          },
          {
            "name": "BM_AddFields/8388608/threads:4",
            "value": 56170078.75000013,
            "unit": "ns/iter",
            "extra": "iterations: 8\ncpu: 113452299.99999993 ns\nthreads: 4"
          },
          {
            "name": "BM_AddFields/8388608/threads:8",
            "value": 60941140.0937514,
            "unit": "ns/iter",
            "extra": "iterations: 8\ncpu: 113586062.50000012 ns\nthreads: 8"
          },
          {
            "name": "BM_MulScalar/8/threads:1",
            "value": 1493.9060254501367,
            "unit": "ns/iter",
            "extra": "iterations: 461061\ncpu: 1493.7635150229592 ns\nthreads: 1"
          },
          {
            "name": "BM_MulScalar/8/threads:2",
            "value": 2160.267160114698,
            "unit": "ns/iter",
            "extra": "iterations: 186916\ncpu: 3742.262299642636 ns\nthreads: 2"
          },
          {
            "name": "BM_MulScalar/8/threads:4",
            "value": 2037.0440127146244,
            "unit": "ns/iter",
            "extra": "iterations: 187656\ncpu: 3717.5038368077803 ns\nthreads: 4"
          },
          {
            "name": "BM_MulScalar/8/threads:8",
            "value": 2101.666785438603,
            "unit": "ns/iter",
            "extra": "iterations: 186632\ncpu: 4165.905632474607 ns\nthreads: 8"
          },
          {
            "name": "BM_MulScalar/64/threads:1",
            "value": 1520.6686951508905,
            "unit": "ns/iter",
            "extra": "iterations: 444174\ncpu: 1520.4856655274737 ns\nthreads: 1"
          },
          {
            "name": "BM_MulScalar/64/threads:2",
            "value": 2178.293633723619,
            "unit": "ns/iter",
            "extra": "iterations: 184582\ncpu: 3839.302315502065 ns\nthreads: 2"
          },
          {
            "name": "BM_MulScalar/64/threads:4",
            "value": 2089.260149010016,
            "unit": "ns/iter",
            "extra": "iterations: 261996\ncpu: 3738.963953648147 ns\nthreads: 4"
          },
          {
            "name": "BM_MulScalar/64/threads:8",
            "value": 2102.8173854755078,
            "unit": "ns/iter",
            "extra": "iterations: 176032\ncpu: 4129.587802217784 ns\nthreads: 8"
          },
          {
            "name": "BM_MulScalar/512/threads:1",
            "value": 1901.9436578527052,
            "unit": "ns/iter",
            "extra": "iterations: 338521\ncpu: 1901.5449558520866 ns\nthreads: 1"
          },
          {
            "name": "BM_MulScalar/512/threads:2",
            "value": 2353.850880141227,
            "unit": "ns/iter",
            "extra": "iterations: 190424\ncpu: 4140.877725496798 ns\nthreads: 2"
          },
          {
            "name": "BM_MulScalar/512/threads:4",
            "value": 2335.0426602396656,
            "unit": "ns/iter",
            "extra": "iterations: 168248\ncpu: 4037.9053540012683 ns\nthreads: 4"
          },
          {
            "name": "BM_MulScalar/512/threads:8",
            "value": 2228.761380963414,
            "unit": "ns/iter",
            "extra": "iterations: 161992\ncpu: 4264.216751444484 ns\nthreads: 8"
          },
          {
            "name": "BM_MulScalar/4096/threads:1",
            "value": 4398.024050600368,
            "unit": "ns/iter",
            "extra": "iterations: 155256\ncpu: 4397.310248879273 ns\nthreads: 1"
          },
          {
            "name": "BM_MulScalar/4096/threads:2",
            "value": 3669.307188430268,
            "unit": "ns/iter",
            "extra": "iterations: 102818\ncpu: 6635.108638565239 ns\nthreads: 2"
          },
          {
            "name": "BM_MulScalar/4096/threads:4",
            "value": 3574.7717787992183,
            "unit": "ns/iter",
            "extra": "iterations: 105768\ncpu: 6356.968081083125 ns\nthreads: 4"
          },
          {
            "name": "BM_MulScalar/4096/threads:8",
            "value": 3326.9828156250014,
            "unit": "ns/iter",
            "extra": "iterations: 80000\ncpu: 6497.212499999938 ns\nthreads: 8"
          },
          {
            "name": "BM_MulScalar/32768/threads:1",
            "value": 27956.657812316284,
            "unit": "ns/iter",
            "extra": "iterations: 25223\ncpu: 27947.952265789052 ns\nthreads: 1"
          },
          {
            "name": "BM_MulScalar/32768/threads:2",
            "value": 14789.021324941332,
            "unit": "ns/iter",
            "extra": "iterations: 23986\ncpu: 28419.14450095902 ns\nthreads: 2"
          },
          {
            "name": "BM_MulScalar/32768/threads:4",
            "value": 13571.562362030962,
            "unit": "ns/iter",
            "extra": "iterations: 25368\ncpu: 27404.679123304937 ns\nthreads: 4"
          },
          {
            "name": "BM_MulScalar/32768/threads:8",
            "value": 14133.128649350496,
            "unit": "ns/iter",
            "extra": "iterations: 22504\ncpu: 30066.89033060778 ns\nthreads: 8"
          },
          {
            "name": "BM_MulScalar/262144/threads:1",
            "value": 284982.40054708044,
            "unit": "ns/iter",
            "extra": "iterations: 2559\ncpu: 284941.2661195773 ns\nthreads: 1"
          },
          {
            "name": "BM_MulScalar/262144/threads:2",
            "value": 138150.0729801986,
            "unit": "ns/iter",
            "extra": "iterations: 2624\ncpu: 275502.5533536606 ns\nthreads: 2"
          },
          {
            "name": "BM_MulScalar/262144/threads:4",
            "value": 137232.9129709599,
            "unit": "ns/iter",
            "extra": "iterations: 2548\ncpu: 277562.2448979572 ns\nthreads: 4"
          },
          {
            "name": "BM_MulScalar/262144/threads:8",
            "value": 135811.82217507713,
            "unit": "ns/iter",
            "extra": "iterations: 2416\ncpu: 287473.8824503323 ns\nthreads: 8"
          },
          {
            "name": "BM_MulScalar/2097152/threads:1",
            "value": 2262782.1743424046,
            "unit": "ns/iter",
            "extra": "iterations: 304\ncpu: 2262002.9605263136 ns\nthreads: 1"
          },
          {
            "name": "BM_MulScalar/2097152/threads:2",
            "value": 1347479.1672931723,
            "unit": "ns/iter",
            "extra": "iterations: 266\ncpu: 2666328.5714285793 ns\nthreads: 2"
          },
          {
            "name": "BM_MulScalar/2097152/threads:4",
            "value": 1698883.4667555925,
            "unit": "ns/iter",
            "extra": "iterations: 188\ncpu: 3481259.0425532004 ns\nthreads: 4"
          },
          {
            "name": "BM_MulScalar/2097152/threads:8",
            "value": 2107582.646875006,
            "unit": "ns/iter",
            "extra": "iterations: 160\ncpu: 4434548.125000004 ns\nthreads: 8"
          },
          {
            "name": "BM_MulScalar/8388608/threads:1",
            "value": 49610502.50000198,
            "unit": "ns/iter",
            "extra": "iterations: 14\ncpu: 49608785.71428548 ns\nthreads: 1"
          },
          {
            "name": "BM_MulScalar/8388608/threads:2",
            "value": 27651557.70833398,
            "unit": "ns/iter",
            "extra": "iterations: 12\ncpu: 54028883.33333335 ns\nthreads: 2"
          },
          {
            "name": "BM_MulScalar/8388608/threads:4",
            "value": 27697414.312501203,
            "unit": "ns/iter",
            "extra": "iterations: 12\ncpu: 56709466.66666695 ns\nthreads: 4"
          },
          {
            "name": "BM_MulScalar/8388608/threads:8",
            "value": 28203596.804689646,
            "unit": "ns/iter",
            "extra": "iterations: 16\ncpu: 58432974.99999997 ns\nthreads: 8"
          },
          {
            "name": "BM_MulFields/8/threads:1",
            "value": 3076.635647602624,
            "unit": "ns/iter",
            "extra": "iterations: 220660\ncpu: 3076.23991661381 ns\nthreads: 1"
          },
          {
            "name": "BM_MulFields/8/threads:2",
            "value": 4787.922323714923,
            "unit": "ns/iter",
            "extra": "iterations: 87472\ncpu: 8330.736692884622 ns\nthreads: 2"
          },
          {
            "name": "BM_MulFields/8/threads:4",
            "value": 4521.406303775293,
            "unit": "ns/iter",
            "extra": "iterations: 111576\ncpu: 8084.297698429799 ns\nthreads: 4"
          },
          {
            "name": "BM_MulFields/8/threads:8",
            "value": 4299.992544345342,
            "unit": "ns/iter",
            "extra": "iterations: 80392\ncpu: 8521.545676186735 ns\nthreads: 8"
          },
          {
            "name": "BM_MulFields/64/threads:1",
            "value": 3120.633646580047,
            "unit": "ns/iter",
            "extra": "iterations: 212931\ncpu: 3119.8209748697755 ns\nthreads: 1"
          },
          {
            "name": "BM_MulFields/64/threads:2",
            "value": 4613.828360944612,
            "unit": "ns/iter",
            "extra": "iterations: 87282\ncpu: 7911.752709607938 ns\nthreads: 2"
          },
          {
            "name": "BM_MulFields/64/threads:4",
            "value": 3963.865744452282,
            "unit": "ns/iter",
            "extra": "iterations: 94808\ncpu: 7038.009450679247 ns\nthreads: 4"
          },
          {
            "name": "BM_MulFields/64/threads:8",
            "value": 4384.233147730101,
            "unit": "ns/iter",
            "extra": "iterations: 89528\ncpu: 8374.965373961226 ns\nthreads: 8"
          },
          {
            "name": "BM_MulFields/512/threads:1",
            "value": 3899.1528421119156,
            "unit": "ns/iter",
            "extra": "iterations: 181344\ncpu: 3897.2510808187985 ns\nthreads: 1"
          },
          {
            "name": "BM_MulFields/512/threads:2",
            "value": 5104.609230921735,
            "unit": "ns/iter",
            "extra": "iterations: 80902\ncpu: 8562.033077056163 ns\nthreads: 2"
          },
          {
            "name": "BM_MulFields/512/threads:4",
            "value": 4857.302818226545,
            "unit": "ns/iter",
            "extra": "iterations: 77052\ncpu: 8401.671598401088 ns\nthreads: 4"
          },
          {
            "name": "BM_MulFields/512/threads:8",
            "value": 4992.893960038898,
            "unit": "ns/iter",
            "extra": "iterations: 92840\ncpu: 9238.675140025805 ns\nthreads: 8"
          },
          {
            "name": "BM_MulFields/4096/threads:1",
            "value": 8675.711133729188,
            "unit": "ns/iter",
            "extra": "iterations: 85479\ncpu: 8631.246271013937 ns\nthreads: 1"
          },
          {
            "name": "BM_MulFields/4096/threads:2",
            "value": 7264.723458188366,
            "unit": "ns/iter",
            "extra": "iterations: 57400\ncpu: 12805.33275261324 ns\nthreads: 2"
          },
          {
            "name": "BM_MulFields/4096/threads:4",
            "value": 7369.325086098461,
            "unit": "ns/iter",
            "extra": "iterations: 54008\ncpu: 13585.900237001928 ns\nthreads: 4"
          },
          {
            "name": "BM_MulFields/4096/threads:8",
            "value": 6484.227670338363,
            "unit": "ns/iter",
            "extra": "iterations: 58560\ncpu: 12688.990778688514 ns\nthreads: 8"
          },
          {
            "name": "BM_MulFields/32768/threads:1",
            "value": 62819.24459692029,
            "unit": "ns/iter",
            "extra": "iterations: 11660\ncpu: 62755.231560892105 ns\nthreads: 1"
          },
          {
            "name": "BM_MulFields/32768/threads:2",
            "value": 32098.647076946876,
            "unit": "ns/iter",
            "extra": "iterations: 11358\ncpu: 63672.90015847849 ns\nthreads: 2"
          },
          {
            "name": "BM_MulFields/32768/threads:4",
            "value": 32050.433748583764,
            "unit": "ns/iter",
            "extra": "iterations: 10596\ncpu: 63605.001887504564 ns\nthreads: 4"
          },
          {
            "name": "BM_MulFields/32768/threads:8",
            "value": 30213.871890625298,
            "unit": "ns/iter",
            "extra": "iterations: 8000\ncpu: 63145.03749999975 ns\nthreads: 8"
          },
          {
            "name": "BM_MulFields/262144/threads:1",
            "value": 482582.7124528278,
            "unit": "ns/iter",
            "extra": "iterations: 1325\ncpu: 482379.3962264127 ns\nthreads: 1"
          },
          {
            "name": "BM_MulFields/262144/threads:2",
            "value": 245837.44123638797,
            "unit": "ns/iter",
            "extra": "iterations: 1472\ncpu: 490921.7391304359 ns\nthreads: 2"
          },
          {
            "name": "BM_MulFields/262144/threads:4",
            "value": 237768.7740778814,
            "unit": "ns/iter",
            "extra": "iterations: 1464\ncpu: 486940.91530054447 ns\nthreads: 4"
          },
          {
            "name": "BM_MulFields/262144/threads:8",
            "value": 263838.021894266,
            "unit": "ns/iter",
            "extra": "iterations: 1296\ncpu: 540447.2222222224 ns\nthreads: 8"
          },
          {
            "name": "BM_MulFields/2097152/threads:1",
            "value": 16718145.374997562,
            "unit": "ns/iter",
            "extra": "iterations: 40\ncpu: 16715542.500000069 ns\nthreads: 1"
          },
          {
            "name": "BM_MulFields/2097152/threads:2",
            "value": 6486896.051723471,
            "unit": "ns/iter",
            "extra": "iterations: 58\ncpu: 12955662.068965571 ns\nthreads: 2"
          },
          {
            "name": "BM_MulFields/2097152/threads:4",
            "value": 5190227.754464364,
            "unit": "ns/iter",
            "extra": "iterations: 56\ncpu: 11924825.000000048 ns\nthreads: 4"
          },
          {
            "name": "BM_MulFields/2097152/threads:8",
            "value": 5089784.587890556,
            "unit": "ns/iter",
            "extra": "iterations: 64\ncpu: 11628196.875 ns\nthreads: 8"
          },
          {
            "name": "BM_MulFields/8388608/threads:1",
            "value": 105638238.14285505,
            "unit": "ns/iter",
            "extra": "iterations: 7\ncpu: 105615528.5714289 ns\nthreads: 1"
          },
          {
            "name": "BM_MulFields/8388608/threads:2",
            "value": 57963603.5000135,
            "unit": "ns/iter",
            "extra": "iterations: 6\ncpu: 114122283.3333324 ns\nthreads: 2"
          },
          {
            "name": "BM_MulFields/8388608/threads:4",
            "value": 57324986.28124816,
            "unit": "ns/iter",
            "extra": "iterations: 8\ncpu: 114944825.00000028 ns\nthreads: 4"
          },
          {
            "name": "BM_MulFields/8388608/threads:8",
            "value": 60294509.81250051,
            "unit": "ns/iter",
            "extra": "iterations: 8\ncpu: 113002587.49999954 ns\nthreads: 8"
          }
        ]
      }
    ]
  }
}
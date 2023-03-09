
// tiltdrivepro_in_neg_table generated by DK/circ_table_gen.py -- do not modify manually

 // variables used
 // --sig_max  -1.500000
 // --table_div  1.700617
 // --table_op  1.000000

struct tabletiltdrivepro_in_neg { // 1-dimensional function table
    float low;
    float high;
    float istep;
    int size;
    float data[];
};

template <int tab_size>
struct tabletiltdrivepro_in_neg_imp {
    float low;
    float high;
    float istep;
    int size;
    float data[tab_size];
    operator tabletiltdrivepro_in_neg&() const { return *(tabletiltdrivepro_in_neg*)this; }
};

 static tabletiltdrivepro_in_neg_imp<2048> tiltdrivepro_in_neg_table __rt_data = {
	0,0.715658,1364.67,2048, {
	0.000000000000,0.001433353226,0.002879693202,0.004337116098,0.005803835880,
	0.007278178271,0.008758574994,0.010243558293,0.011731755709,0.013221885110,
	0.014712749956,0.016203234788,0.017692300938,0.019178982442,0.020662382154,
	0.022141668041,0.023616069663,0.025084874822,0.026547426373,0.028003119195,
	0.029451397302,0.030891751105,0.032323714802,0.033746863902,0.035160812865,
	0.036565212868,0.037959749672,0.039344141606,0.040718137644,0.042081515585,
	0.043434080319,0.044775662185,0.046106115412,0.047425316641,0.048733163516,
	0.050029573357,0.051314481898,0.052587842084,0.053849622944,0.055099808511,
	0.056338396806,0.057565398869,0.058780837862,0.059984748179,0.061177174650,
	0.062358171758,0.063527802908,0.064686139742,0.065833261479,0.066969254300,
	0.068094210770,0.069208229282,0.070311413544,0.071403872086,0.072485717803,
	0.073557067520,0.074618041581,0.075668763469,0.076709359440,0.077739958188,
	0.078760690520,0.079771689064,0.080773087981,0.081765022707,0.082747629703,
	0.083721046227,0.084685410117,0.085640859587,0.086587533047,0.087525568918,
	0.088455105478,0.089376280703,0.090289232131,0.091194096732,0.092091010784,
	0.092980109764,0.093861528244,0.094735399797,0.095601856910,0.096461030904,
	0.097313051864,0.098158048568,0.098996148431,0.099827477447,0.100652160147,
	0.101470319546,0.102282077113,0.103087552729,0.103886864664,0.104680129544,
	0.105467462335,0.106248976321,0.107024783087,0.107794992510,0.108559712747,
	0.109319050230,0.110073109660,0.110821994006,0.111565804505,0.112304640664,
	0.113038600265,0.113767779369,0.114492272327,0.115212171784,0.115927568696,
	0.116638552334,0.117345210305,0.118047628558,0.118745891407,0.119440081539,
	0.120130280037,0.120816566393,0.121499018530,0.122177712819,0.122852724094,
	0.123524125682,0.124191989411,0.124856385641,0.125517383277,0.126175049796,
	0.126829451264,0.127480652360,0.128128716397,0.128773705344,0.129415679848,
	0.130054699256,0.130690821634,0.131324103794,0.131954601313,0.132582368556,
	0.133207458694,0.133829923733,0.134449814528,0.135067180809,0.135682071200,
	0.136294533243,0.136904613413,0.137512357147,0.138117808858,0.138721011956,
	0.139322008873,0.139920841075,0.140517549089,0.141112172519,0.141704750063,
	0.142295319536,0.142883917887,0.143470581216,0.144055344793,0.144638243076,
	0.145219309729,0.145798577637,0.146376078925,0.146951844976,0.147525906442,
	0.148098293266,0.148669034695,0.149238159295,0.149805694967,0.150371668964,
	0.150936107903,0.151499037778,0.152060483980,0.152620471305,0.153179023971,
	0.153736165629,0.154291919379,0.154846307781,0.155399352867,0.155951076154,
	0.156501498657,0.157050640900,0.157598522927,0.158145164316,0.158690584184,
	0.159234801205,0.159777833618,0.160319699235,0.160860415454,0.161399999268,
	0.161938467275,0.162475835687,0.163012120341,0.163547336704,0.164081499887,
	0.164614624651,0.165146725415,0.165677816267,0.166207910967,0.166737022962,
	0.167265165389,0.167792351081,0.168318592579,0.168843902139,0.169368291734,
	0.169891773065,0.170414357569,0.170936056421,0.171456880544,0.171976840614,
	0.172495947065,0.173014210100,0.173531639688,0.174048245579,0.174564037302,
	0.175079024176,0.175593215312,0.176106619619,0.176619245810,0.177131102406,
	0.177642197739,0.178152539961,0.178662137045,0.179170996790,0.179679126826,
	0.180186534619,0.180693227474,0.181199212538,0.181704496806,0.182209087124,
	0.182712990193,0.183216212572,0.183718760680,0.184220640805,0.184721859100,
	0.185222421592,0.185722334182,0.186221602651,0.186720232660,0.187218229753,
	0.187715599364,0.188212346815,0.188708477320,0.189203995990,0.189698907833,
	0.190193217758,0.190686930577,0.191180051007,0.191672583672,0.192164533107,
	0.192655903759,0.193146699989,0.193636926075,0.194126586213,0.194615684519,
	0.195104225035,0.195592211722,0.196079648471,0.196566539099,0.197052887356,
	0.197538696920,0.198023971404,0.198508714355,0.198992929257,0.199476619533,
	0.199959788543,0.200442439589,0.200924575917,0.201406200715,0.201887317116,
	0.202367928200,0.202848036995,0.203327646478,0.203806759575,0.204285379165,
	0.204763508077,0.205241149098,0.205718304965,0.206194978373,0.206671171975,
	0.207146888380,0.207622130157,0.208096899834,0.208571199900,0.209045032806,
	0.209518400964,0.209991306753,0.210463752513,0.210935740550,0.211407273137,
	0.211878352511,0.212348980879,0.212819160415,0.213288893262,0.213758181534,
	0.214227027312,0.214695432652,0.215163399578,0.215630930089,0.216098026154,
	0.216564689719,0.217030922700,0.217496726991,0.217962104459,0.218427056947,
	0.218891586276,0.219355694241,0.219819382617,0.220282653153,0.220745507581,
	0.221207947607,0.221669974918,0.222131591182,0.222592798045,0.223053597134,
	0.223513990057,0.223973978402,0.224433563740,0.224892747623,0.225351531586,
	0.225809917146,0.226267905803,0.226725499042,0.227182698329,0.227639505115,
	0.228095920836,0.228551946913,0.229007584750,0.229462835737,0.229917701251,
	0.230372182652,0.230826281288,0.231279998493,0.231733335588,0.232186293879,
	0.232638874660,0.233091079213,0.233542908808,0.233994364699,0.234445448132,
	0.234896160341,0.235346502545,0.235796475954,0.236246081768,0.236695321173,
	0.237144195346,0.237592705453,0.238040852650,0.238488638081,0.238936062882,
	0.239383128178,0.239829835084,0.240276184706,0.240722178139,0.241167816472,
	0.241613100781,0.242058032136,0.242502611595,0.242946840211,0.243390719025,
	0.243834249072,0.244277431377,0.244720266958,0.245162756823,0.245604901974,
	0.246046703406,0.246488162102,0.246929279042,0.247370055196,0.247810491527,
	0.248250588991,0.248690348536,0.249129771105,0.249568857631,0.250007609043,
	0.250446026260,0.250884110198,0.251321861763,0.251759281857,0.252196371373,
	0.252633131200,0.253069562220,0.253505665308,0.253941441333,0.254376891160,
	0.254812015644,0.255246815639,0.255681291988,0.256115445533,0.256549277107,
	0.256982787538,0.257415977651,0.257848848261,0.258281400182,0.258713634220,
	0.259145551176,0.259577151847,0.260008437023,0.260439407490,0.260870064029,
	0.261300407415,0.261730438420,0.262160157809,0.262589566344,0.263018664779,
	0.263447453867,0.263875934355,0.264304106984,0.264731972493,0.265159531613,
	0.265586785073,0.266013733598,0.266440377906,0.266866718713,0.267292756730,
	0.267718492663,0.268143927214,0.268569061081,0.268993894959,0.269418429536,
	0.269842665498,0.270266603528,0.270690244302,0.271113588494,0.271536636773,
	0.271959389806,0.272381848254,0.272804012775,0.273225884023,0.273647462649,
	0.274068749299,0.274489744616,0.274910449239,0.275330863805,0.275750988944,
	0.276170825285,0.276590373454,0.277009634071,0.277428607755,0.277847295119,
	0.278265696774,0.278683813329,0.279101645386,0.279519193547,0.279936458410,
	0.280353440568,0.280770140612,0.281186559130,0.281602696705,0.282018553920,
	0.282434131352,0.282849429576,0.283264449164,0.283679190684,0.284093654701,
	0.284507841778,0.284921752474,0.285335387347,0.285748746948,0.286161831828,
	0.286574642536,0.286987179614,0.287399443606,0.287811435049,0.288223154479,
	0.288634602430,0.289045779432,0.289456686011,0.289867322692,0.290277689998,
	0.290687788447,0.291097618555,0.291507180836,0.291916475802,0.292325503959,
	0.292734265814,0.293142761870,0.293550992626,0.293958958582,0.294366660230,
	0.294774098065,0.295181272576,0.295588184251,0.295994833574,0.296401221029,
	0.296807347094,0.297213212247,0.297618816963,0.298024161716,0.298429246974,
	0.298834073207,0.299238640878,0.299642950451,0.300047002387,0.300450797144,
	0.300854335178,0.301257616942,0.301660642887,0.302063413464,0.302465929118,
	0.302868190293,0.303270197433,0.303671950977,0.304073451363,0.304474699027,
	0.304875694401,0.305276437918,0.305676930006,0.306077171091,0.306477161600,
	0.306876901954,0.307276392574,0.307675633878,0.308074626283,0.308473370203,
	0.308871866051,0.309270114235,0.309668115165,0.310065869247,0.310463376884,
	0.310860638478,0.311257654431,0.311654425139,0.312050951000,0.312447232406,
	0.312843269752,0.313239063426,0.313634613817,0.314029921312,0.314424986295,
	0.314819809149,0.315214390255,0.315608729992,0.316002828737,0.316396686865,
	0.316790304749,0.317183682763,0.317576821274,0.317969720651,0.318362381261,
	0.318754803468,0.319146987634,0.319538934121,0.319930643287,0.320322115491,
	0.320713351087,0.321104350431,0.321495113873,0.321885641765,0.322275934456,
	0.322665992293,0.323055815621,0.323445404784,0.323834760125,0.324223881985,
	0.324612770701,0.325001426612,0.325389850053,0.325778041358,0.326166000861,
	0.326553728891,0.326941225779,0.327328491852,0.327715527436,0.328102332856,
	0.328488908435,0.328875254496,0.329261371357,0.329647259337,0.330032918755,
	0.330418349924,0.330803553160,0.331188528775,0.331573277080,0.331957798384,
	0.332342092997,0.332726161224,0.333110003371,0.333493619743,0.333877010640,
	0.334260176365,0.334643117218,0.335025833495,0.335408325495,0.335790593512,
	0.336172637842,0.336554458775,0.336936056605,0.337317431620,0.337698584110,
	0.338079514362,0.338460222661,0.338840709293,0.339220974540,0.339601018686,
	0.339980842010,0.340360444791,0.340739827309,0.341118989840,0.341497932659,
	0.341876656042,0.342255160260,0.342633445585,0.343011512289,0.343389360641,
	0.343766990908,0.344144403358,0.344521598256,0.344898575866,0.345275336452,
	0.345651880276,0.346028207598,0.346404318679,0.346780213776,0.347155893148,
	0.347531357049,0.347906605735,0.348281639461,0.348656458477,0.349031063037,
	0.349405453390,0.349779629786,0.350153592472,0.350527341697,0.350900877704,
	0.351274200741,0.351647311049,0.352020208872,0.352392894452,0.352765368028,
	0.353137629840,0.353509680127,0.353881519125,0.354253147070,0.354624564199,
	0.354995770744,0.355366766940,0.355737553017,0.356108129207,0.356478495740,
	0.356848652845,0.357218600749,0.357588339680,0.357957869864,0.358327191525,
	0.358696304887,0.359065210173,0.359433907606,0.359802397406,0.360170679793,
	0.360538754986,0.360906623204,0.361274284663,0.361641739580,0.362008988169,
	0.362376030646,0.362742867223,0.363109498114,0.363475923529,0.363842143679,
	0.364208158773,0.364573969022,0.364939574632,0.365304975810,0.365670172763,
	0.366035165696,0.366399954813,0.366764540317,0.367128922412,0.367493101298,
	0.367857077177,0.368220850248,0.368584420711,0.368947788764,0.369310954604,
	0.369673918428,0.370036680432,0.370399240811,0.370761599758,0.371123757467,
	0.371485714130,0.371847469940,0.372209025086,0.372570379758,0.372931534147,
	0.373292488440,0.373653242825,0.374013797489,0.374374152618,0.374734308396,
	0.375094265010,0.375454022641,0.375813581474,0.376172941690,0.376532103471,
	0.376891066997,0.377249832449,0.377608400005,0.377966769844,0.378324942143,
	0.378682917080,0.379040694831,0.379398275571,0.379755659475,0.380112846717,
	0.380469837470,0.380826631907,0.381183230200,0.381539632520,0.381895839037,
	0.382251849922,0.382607665343,0.382963285468,0.383318710467,0.383673940505,
	0.384028975748,0.384383816364,0.384738462516,0.385092914369,0.385447172086,
	0.385801235831,0.386155105767,0.386508782053,0.386862264853,0.387215554326,
	0.387568650632,0.387921553930,0.388274264378,0.388626782135,0.388979107358,
	0.389331240203,0.389683180826,0.390034929382,0.390386486027,0.390737850915,
	0.391089024198,0.391440006030,0.391790796563,0.392141395949,0.392491804340,
	0.392842021884,0.393192048733,0.393541885036,0.393891530941,0.394240986597,
	0.394590252151,0.394939327750,0.395288213541,0.395636909670,0.395985416281,
	0.396333733520,0.396681861530,0.397029800456,0.397377550441,0.397725111626,
	0.398072484154,0.398419668167,0.398766663804,0.399113471207,0.399460090516,
	0.399806521869,0.400152765405,0.400498821263,0.400844689579,0.401190370493,
	0.401535864139,0.401881170654,0.402226290174,0.402571222833,0.402915968767,
	0.403260528109,0.403604900993,0.403949087552,0.404293087918,0.404636902224,
	0.404980530601,0.405323973181,0.405667230093,0.406010301468,0.406353187436,
	0.406695888125,0.407038403664,0.407380734182,0.407722879807,0.408064840664,
	0.408406616882,0.408748208586,0.409089615903,0.409430838957,0.409771877874,
	0.410112732778,0.410453403793,0.410793891043,0.411134194650,0.411474314738,
	0.411814251428,0.412154004843,0.412493575103,0.412832962330,0.413172166644,
	0.413511188165,0.413850027012,0.414188683305,0.414527157163,0.414865448704,
	0.415203558045,0.415541485304,0.415879230598,0.416216794045,0.416554175759,
	0.416891375857,0.417228394454,0.417565231665,0.417901887606,0.418238362389,
	0.418574656129,0.418910768938,0.419246700931,0.419582452220,0.419918022917,
	0.420253413133,0.420588622980,0.420923652569,0.421258502011,0.421593171416,
	0.421927660894,0.422261970554,0.422596100505,0.422930050857,0.423263821717,
	0.423597413194,0.423930825395,0.424264058428,0.424597112398,0.424929987414,
	0.425262683581,0.425595201004,0.425927539790,0.426259700044,0.426591681869,
	0.426923485371,0.427255110654,0.427586557820,0.427917826975,0.428248918220,
	0.428579831658,0.428910567392,0.429241125523,0.429571506153,0.429901709384,
	0.430231735316,0.430561584051,0.430891255687,0.431220750326,0.431550068067,
	0.431879209009,0.432208173251,0.432536960892,0.432865572030,0.433194006763,
	0.433522265189,0.433850347406,0.434178253509,0.434505983597,0.434833537766,
	0.435160916111,0.435488118728,0.435815145714,0.436141997162,0.436468673169,
	0.436795173828,0.437121499235,0.437447649482,0.437773624664,0.438099424874,
	0.438425050205,0.438750500750,0.439075776602,0.439400877852,0.439725804593,
	0.440050556917,0.440375134914,0.440699538675,0.441023768293,0.441347823856,
	0.441671705456,0.441995413182,0.442318947124,0.442642307371,0.442965494012,
	0.443288507137,0.443611346834,0.443934013191,0.444256506296,0.444578826237,
	0.444900973101,0.445222946976,0.445544747949,0.445866376106,0.446187831534,
	0.446509114319,0.446830224547,0.447151162303,0.447471927673,0.447792520742,
	0.448112941596,0.448433190318,0.448753266993,0.449073171705,0.449392904538,
	0.449712465576,0.450031854903,0.450351072600,0.450670118751,0.450988993440,
	0.451307696747,0.451626228756,0.451944589549,0.452262779206,0.452580797810,
	0.452898645442,0.453216322182,0.453533828112,0.453851163312,0.454168327862,
	0.454485321843,0.454802145334,0.455118798414,0.455435281164,0.455751593662,
	0.456067735988,0.456383708219,0.456699510435,0.457015142713,0.457330605133,
	0.457645897771,0.457961020704,0.458275974012,0.458590757770,0.458905372056,
	0.459219816946,0.459534092517,0.459848198846,0.460162136007,0.460475904078,
	0.460789503134,0.461102933250,0.461416194502,0.461729286964,0.462042210712,
	0.462354965820,0.462667552362,0.462979970414,0.463292220048,0.463604301339,
	0.463916214360,0.464227959185,0.464539535887,0.464850944540,0.465162185215,
	0.465473257985,0.465784162924,0.466094900103,0.466405469595,0.466715871471,
	0.467026105802,0.467336172661,0.467646072119,0.467955804247,0.468265369116,
	0.468574766796,0.468883997359,0.469193060874,0.469501957412,0.469810687043,
	0.470119249836,0.470427645862,0.470735875189,0.471043937888,0.471351834026,
	0.471659563674,0.471967126900,0.472274523772,0.472581754359,0.472888818730,
	0.473195716951,0.473502449092,0.473809015220,0.474115415402,0.474421649706,
	0.474727718199,0.475033620948,0.475339358020,0.475644929481,0.475950335399,
	0.476255575839,0.476560650869,0.476865560553,0.477170304958,0.477474884149,
	0.477779298192,0.478083547153,0.478387631097,0.478691550088,0.478995304192,
	0.479298893474,0.479602317998,0.479905577828,0.480208673029,0.480511603665,
	0.480814369800,0.481116971498,0.481419408822,0.481721681837,0.482023790604,
	0.482325735189,0.482627515653,0.482929132060,0.483230584472,0.483531872952,
	0.483832997563,0.484133958367,0.484434755426,0.484735388802,0.485035858557,
	0.485336164752,0.485636307450,0.485936286712,0.486236102599,0.486535755173,
	0.486835244493,0.487134570622,0.487433733620,0.487732733548,0.488031570466,
	0.488330244435,0.488628755514,0.488927103764,0.489225289245,0.489523312016,
	0.489821172137,0.490118869668,0.490416404668,0.490713777197,0.491010987313,
	0.491308035076,0.491604920544,0.491901643776,0.492198204831,0.492494603766,
	0.492790840642,0.493086915515,0.493382828444,0.493678579487,0.493974168702,
	0.494269596145,0.494564861876,0.494859965951,0.495154908428,0.495449689364,
	0.495744308815,0.496038766840,0.496333063494,0.496627198835,0.496921172919,
	0.497214985803,0.497508637542,0.497802128194,0.498095457814,0.498388626458,
	0.498681634182,0.498974481042,0.499267167094,0.499559692393,0.499852056994,
	0.500144260953,0.500436304324,0.500728187164,0.501019909526,0.501311471467,
	0.501602873039,0.501894114299,0.502185195300,0.502476116097,0.502766876744,
	0.503057477296,0.503347917806,0.503638198328,0.503928318916,0.504218279624,
	0.504508080506,0.504797721614,0.505087203002,0.505376524724,0.505665686833,
	0.505954689382,0.506243532423,0.506532216010,0.506820740195,0.507109105031,
	0.507397310571,0.507685356866,0.507973243970,0.508260971934,0.508548540811,
	0.508835950653,0.509123201510,0.509410293437,0.509697226483,0.509984000701,
	0.510270616142,0.510557072857,0.510843370899,0.511129510317,0.511415491164,
	0.511701313489,0.511986977345,0.512272482782,0.512557829850,0.512843018601,
	0.513128049084,0.513412921351,0.513697635451,0.513982191435,0.514266589353,
	0.514550829254,0.514834911190,0.515118835209,0.515402601361,0.515686209697,
	0.515969660266,0.516252953117,0.516536088299,0.516819065862,0.517101885856,
	0.517384548328,0.517667053329,0.517949400907,0.518231591111,0.518513623990,
	0.518795499592,0.519077217967,0.519358779161,0.519640183225,0.519921430206,
	0.520202520152,0.520483453112,0.520764229133,0.521044848264,0.521325310552,
	0.521605616046,0.521885764792,0.522165756839,0.522445592234,0.522725271025,
	0.523004793258,0.523284158982,0.523563368243,0.523842421088,0.524121317564,
	0.524400057720,0.524678641600,0.524957069252,0.525235340723,0.525513456059,
	0.525791415308,0.526069218514,0.526346865725,0.526624356987,0.526901692346,
	0.527178871849,0.527455895541,0.527732763468,0.528009475676,0.528286032211,
	0.528562433119,0.528838678445,0.529114768236,0.529390702535,0.529666481390,
	0.529942104844,0.530217572944,0.530492885735,0.530768043261,0.531043045569,
	0.531317892701,0.531592584705,0.531867121624,0.532141503503,0.532415730386,
	0.532689802319,0.532963719346,0.533237481512,0.533511088860,0.533784541434,
	0.534057839280,0.534330982441,0.534603970962,0.534876804886,0.535149484257,
	0.535422009119,0.535694379516,0.535966595492,0.536238657090,0.536510564353,
	0.536782317326,0.537053916051,0.537325360573,0.537596650933,0.537867787177,
	0.538138769345,0.538409597483,0.538680271632,0.538950791836,0.539221158137,
	0.539491370578,0.539761429203,0.540031334053,0.540301085171,0.540570682600,
	0.540840126383,0.541109416560,0.541378553176,0.541647536272,0.541916365890,
	0.542185042072,0.542453564861,0.542721934299,0.542990150427,0.543258213287,
	0.543526122921,0.543793879370,0.544061482678,0.544328932884,0.544596230030,
	0.544863374159,0.545130365311,0.545397203527,0.545663888850,0.545930421319,
	0.546196800977,0.546463027864,0.546729102022,0.546995023491,0.547260792312,
	0.547526408527,0.547791872175,0.548057183297,0.548322341935,0.548587348128,
	0.548852201918,0.549116903344,0.549381452448,0.549645849269,0.549910093847,
	0.550174186224,0.550438126439,0.550701914531,0.550965550542,0.551229034511,
	0.551492366478,0.551755546483,0.552018574566,0.552281450766,0.552544175123,
	0.552806747677,0.553069168467,0.553331437533,0.553593554914,0.553855520650,
	0.554117334780,0.554378997342,0.554640508377,0.554901867924,0.555163076021,
	0.555424132708,0.555685038023,0.555945792006,0.556206394696,0.556466846130,
	0.556727146349,0.556987295390,0.557247293292,0.557507140094,0.557766835835,
	0.558026380552,0.558285774285,0.558545017071,0.558804108949,0.559063049957,
	0.559321840134,0.559580479517,0.559838968144,0.560097306055,0.560355493285,
	0.560613529875,0.560871415860,0.561129151280,0.561386736171,0.561644170573,
	0.561901454521,0.562158588054,0.562415571210,0.562672404026,0.562929086538,
	0.563185618786,0.563442000805,0.563698232634,0.563954314310,0.564210245869,
	0.564466027349,0.564721658786,0.564977140219,0.565232471684,0.565487653217,
	0.565742684856,0.565997566637,0.566252298598,0.566506880774,0.566761313203,
	0.567015595921,0.567269728965,0.567523712372,0.567777546176,0.568031230416,
	0.568284765128,0.568538150347,0.568791386110,0.569044472453,0.569297409412,
	0.569550197024,0.569802835325,0.570055324350,0.570307664135,0.570559854717,
	0.570811896131,0.571063788413,0.571315531598,0.571567125723,0.571818570824,
	0.572069866935,0.572321014092,0.572572012331,0.572822861687,0.573073562197,
	0.573324113894,0.573574516814,0.573824770993,0.574074876466,0.574324833268,
	0.574574641435,0.574824301000,0.575073812000,0.575323174469,0.575572388442,
	0.575821453955,0.576070371041,0.576319139737,0.576567760075,0.576816232092,
	0.577064555822,0.577312731299,0.577560758559,0.577808637635,0.578056368562,
	0.578303951374,0.578551386107,0.578798672793,0.579045811468,0.579292802166,
	0.579539644921,0.579786339767,0.580032886738,0.580279285869,0.580525537192,
	0.580771640743,0.581017596555,0.581263404663,0.581509065099,0.581754577898,
	0.581999943093,0.582245160719,0.582490230808,0.582735153395,0.582979928512,
	0.583224556195,0.583469036475,0.583713369387,0.583957554964,0.584201593238,
	0.584445484245,0.584689228016,0.584932824585,0.585176273985,0.585419576250,
	0.585662731412,0.585905739505,0.586148600561,0.586391314613,0.586633881695,
	0.586876301839,0.587118575078,0.587360701444,0.587602680972,0.587844513692,
	0.588086199639,0.588327738844,0.588569131340,0.588810377160,0.589051476336,
	0.589292428900,0.589533234886,0.589773894325,0.590014407249,0.590254773692,
	0.590494993685,0.590735067260,0.590974994450,0.591214775286,0.591454409801,
	0.591693898027,0.591933239996,0.592172435739,0.592411485289,0.592650388678,
	0.592889145937,0.593127757098,0.593366222192,0.593604541253,0.593842714310,
	0.594080741396,0.594318622543,0.594556357781,0.594793947143,0.595031390660,
	0.595268688363,0.595505840284,0.595742846454,0.595979706904,0.596216421666,
	0.596452990770,0.596689414249,0.596925692132,0.597161824452,0.597397811239,
	0.597633652524,0.597869348338,0.598104898713,0.598340303678,0.598575563265,
	0.598810677506,0.599045646429,0.599280470067,0.599515148450,0.599749681608,
	0.599984069573,0.600218312374,0.600452410043,0.600686362609,0.600920170104,
	0.601153832558,0.601387350000,0.601620722463,0.601853949974,0.602087032566,
	0.602319970268,0.602552763111,0.602785411124,0.603017914338,0.603250272783,
	0.603482486488,0.603714555484,0.603946479801,0.604178259469,0.604409894518,
	0.604641384976,0.604872730875,0.605103932245,0.605334989113,0.605565901512,
	0.605796669469,0.606027293015,0.606257772179,0.606488106991,0.606718297481,
	0.606948343677,0.607178245610,0.607408003308,0.607637616801,0.607867086118,
	0.608096411289,0.608325592343,0.608554629308,0.608783522215,0.609012271093,
	0.609240875969,0.609469336874,0.609697653837,0.609925826886,0.610153856050,
	0.610381741359,0.610609482841,0.610837080525,0.611064534439,0.611291844614,
	0.611519011076,0.611746033856,0.611972912981,0.612199648480,0.612426240382,
	0.612652688716,0.612878993509,0.613105154791,0.613331172589,0.613557046933,
	0.613782777850,0.614008365368,0.614233809517,0.614459110325,0.614684267818,
	0.614909282027,0.615134152978,0.615358880701,0.615583465222,0.615807906571,
	0.616032204774,0.616256359861,0.616480371858,0.616704240795,0.616927966698,
	0.617151549595,0.617374989515,0.617598286485,0.617821440532,0.618044451685,
	0.618267319971,0.618490045418,0.618712628053,0.618935067904,0.619157364998,
	0.619379519362,0.619601531025,0.619823400013,0.620045126354,0.620266710076,
	0.620488151204,0.620709449768,0.620930605793,0.621151619308,0.621372490339,
	0.621593218913,0.621813805057,0.622034248799,0.622254550166,0.622474709184,
	0.622694725880,0.622914600282,0.623134332415,0.623353922308,0.623573369986,
	0.623792675477,0.624011838808,0.624230860004,0.624449739092,0.624668476100,
	0.624887071054,0.625105523980,0.625323834905,0.625542003855,0.625760030857,
	0.625977915937,0.626195659121,0.626413260437,0.626630719910,0.626848037566,
	0.627065213432,0.627282247534,0.627499139898,0.627715890550,0.627932499517,
	0.628148966824,0.628365292498,0.628581476563,0.628797519048,0.629013419977,
	0.629229179376,0.629444797271,0.629660273688,0.629875608652,0.630090802191,
	0.630305854328,0.630520765091,0.630735534503,0.630950162593,0.631164649384,
	0.631378994902,0.631593199173,0.631807262223,0.632021184076,0.632234964759,
	0.632448604296,0.632662102714,0.632875460036,0.633088676290,0.633301751499,
	0.633514685689,0.633727478886,0.633940131114,0.634152642398,0.634365012764,
	0.634577242237,0.634789330842,0.635001278603,0.635213085546,0.635424751695,
	0.635636277076,0.635847661714,0.636058905632,0.636270008857,0.636480971412,
	0.636691793322,0.636902474613,0.637113015308,0.637323415433,0.637533675012,
	0.637743794069,0.637953772629,0.638163610717,0.638373308356,0.638582865573,
	0.638792282389,0.639001558831,0.639210694923,0.639419690688,0.639628546152,
	0.639837261337,0.640045836269,0.640254270972,0.640462565469,0.640670719785,
	0.640878733944,0.641086607970,0.641294341887,0.641501935718,0.641709389489,
	0.641916703222,0.642123876942,0.642330910672,0.642537804436,0.642744558258,
	0.642951172162,0.643157646171,0.643363980309,0.643570174600,0.643776229067,
	0.643982143734,0.644187918624,0.644393553762,0.644599049169,0.644804404871,
	0.645009620889,0.645214697248,0.645419633972,0.645624431082,0.645829088603,
	0.646033606557,0.646237984969,0.646442223861,0.646646323255,0.646850283177,
	0.647054103648,0.647257784691,0.647461326330,0.647664728588,0.647867991486,
	0.648071115050,0.648274099301,0.648476944261,0.648679649955,0.648882216405,
	0.649084643633,0.649286931662,0.649489080516,0.649691090216,0.649892960785,
	0.650094692246,0.650296284621,0.650497737934,0.650699052206,0.650900227459,
	0.651101263717,0.651302161002,0.651502919336,0.651703538742,0.651904019241,
	0.652104360857,0.652304563610,0.652504627525,0.652704552622,0.652904338923,
	0.653103986452,0.653303495230,0.653502865279,0.653702096621,0.653901189278,
	0.654100143273,0.654298958626,0.654497635361,0.654696173498,0.654894573059,
	0.655092834067,0.655290956544,0.655488940510,0.655686785988,0.655884492999,
	0.656082061565,0.656279491707,0.656476783448,0.656673936809,0.656870951810,
	0.657067828474,0.657264566823,0.657461166876,0.657657628657,0.657853952186,
	0.658050137485,0.658246184574,0.658442093476,0.658637864211,0.658833496800,
	0.659028991266,0.659224347628,0.659419565908,0.659614646127,0.659809588306,
	0.660004392467,0.660199058629,0.660393586814,0.660587977044,0.660782229338,
	0.660976343718,0.661170320204,0.661364158817,0.661557859579,0.661751422509,
	0.661944847629,0.662138134958,0.662331284519,0.662524296330,0.662717170414,
	0.662909906790,0.663102505479,0.663294966501,0.663487289878,0.663679475628,
	0.663871523773,0.664063434333,0.664255207328,0.664446842779,0.664638340705,
	0.664829701128,0.665020924067,0.665212009542,0.665402957574,0.665593768182,
	0.665784441387,0.665974977209,0.666165375667,0.666355636782,0.666545760574,
	0.666735747062,0.666925596267,0.667115308208,0.667304882905,0.667494320377,
	0.667683620646,0.667872783729,0.668061809648,0.668250698421,0.668439450069,
	0.668628064610,0.668816542065,0.669004882453,0.669193085793,0.669381152105,
	0.669569081409,0.669756873724,0.669944529068,0.670132047462,0.670319428925,
	0.670506673476,0.670693781135,0.670880751920,0.671067585851,0.671254282946,
	0.671440843226,0.671627266709,0.671813553414,0.671999703361,0.672185716567,
	0.672371593053,0.672557332838,0.672742935939,0.672928402376,0.673113732168,
	0.673298925334,0.673483981892,0.673668901861,0.673853685260,0.674038332108,
	0.674222842423,0.674407216224,0.674591453529,0.674775554358,0.674959518728,
	0.675143346658,0.675327038167,0.675510593272,0.675694011994,0.675877294349,
	0.676060440356,0.676243450033,0.676426323400,0.676609060473,0.676791661272,
	0.676974125814,0.677156454118,0.677338646202,0.677520702083,0.677702621780,
	0.677884405312,0.678066052695,0.678247563949,0.678428939090,0.678610178137,
	0.678791281108,0.678972248020,0.679153078892,0.679333773741,0.679514332585,
	0.679694755442,0.679875042329,0.680055193265,0.680235208266,0.680415087350,
	0.680594830535,0.680774437839,0.680953909279,0.681133244872,0.681312444636,
	0.681491508588,0.681670436746,0.681849229127,0.682027885749,0.682206406628,
	0.682384791783,0.682563041229,0.682741154985,0.682919133068,0.683096975495,
	0.683274682282,0.683452253447,0.683629689008,0.683806988980,0.683984153382,
	0.684161182229,0.684338075539,0.684514833330,0.684691455617,0.684867942417,
	0.685044293748,0.685220509626,0.685396590068,0.685572535090,0.685748344710,
	0.685924018944,0.686099557808,0.686274961319,0.686450229494,0.686625362350,
	0.686800359902,0.686975222167,0.687149949162,0.687324540903,0.687498997406,
	0.687673318689,0.687847504766,0.688021555655,0.688195471371,0.688369251931,
	0.688542897351,0.688716407648,0.688889782837,0.689063022934,0.689236127956,
	0.689409097918,0.689581932837,0.689754632728,0.689927197608,0.690099627493,
	0.690271922397,0.690444082338,0.690616107331,0.690787997391,0.690959752535,
	0.691131372778,0.691302858137,0.691474208626,0.691645424261,0.691816505058,
	0.691987451032,0.692158262200,0.692328938576,0.692499480176,0.692669887015,
	0.692840159109,0.693010296474,0.693180299124,0.693350167075,0.693519900343,
	0.693689498942,0.693858962888,0.694028292196,0.694197486881,0.694366546958,
	0.694535472443,0.694704263351,0.694872919695,0.695041441493,0.695209828758,
	0.695378081506,0.695546199751,0.695714183509,0.695882032794,0.696049747621,
	0.696217328005,0.696384773961,0.696552085504,0.696719262647,0.696886305406,
	0.697053213796,0.697219987831,0.697386627526,0.697553132896,0.697719503954,
	0.697885740715,0.698051843194,0.698217811406,0.698383645364,0.698549345084,
	0.698714910579,0.698880341863,0.699045638952,0.699210801859,0.699375830598,
	0.699540725185,0.699705485632,0.699870111954,0.700034604166,0.700198962280,
	0.700363186312,0.700527276275,0.700691232184,0.700855054051,0.701018741892,
	0.701182295720,0.701345715549,0.701509001393,0.701672153265,0.701835171179,
	0.701998055150,0.702160805190,0.702323421314,0.702485903534,0.702648251866,
	0.702810466322,0.702972546915,0.703134493660,0.703296306570,0.703457985658,
	0.703619530938,0.703780942423,0.703942220127,0.704103364063,0.704264374243,
	0.704425250682,0.704585993393,0.704746602389,0.704907077683,0.705067419289,
	0.705227627218,0.705387701486,0.705547642103,0.705707449085,0.705867122443,
	0.706026662191,0.706186068341,0.706345340906,0.706504479900,0.706663485336,
	0.706822357225,0.706981095581,0.707139700416,0.707298171744,0.707456509577,
	0.707614713928,0.707772784809,0.707930722232,0.708088526212,0.708246196759,
	0.708403733887,0.708561137608,0.708718407934,0.708875544878,0.709032548453,
	0.709189418670,0.709346155542,0.709502759082,0.709659229301,0.709815566212,
	0.709971769827,0.710127840158,0.710283777217,0.710439581017,0.710595251570,
	0.710750788887,0.710906192981,0.711061463863,0.711216601547,0.711371606042,
	0.711526477362,0.711681215519,0.711835820523,0.711990292388,0.712144631124,
	0.712298836744,0.712452909259,0.712606848681,0.712760655022,0.712914328292,
	0.713067868505,0.713221275671,0.713374549802,0.713527690909,0.713680699005,
	0.713833574099,0.713986316205,0.714138925332,0.714291401493,0.714443744699,
	0.714595954961,0.714748032290,0.714899976698,0.715051788196,0.715203466795,
	0.715355012506,0.715506425340,0.715657705309
	}
};

double always_inline tiltdrivepro_in_negclip(double x) {
    double f = fabs(x);
    f = f * tiltdrivepro_in_neg_table.istep;
    int i = static_cast<int>(f);
    if (i < 0) {
        f = tiltdrivepro_in_neg_table.data[0];
    } else if (i >= tiltdrivepro_in_neg_table.size-1) {
        f = tiltdrivepro_in_neg_table.data[tiltdrivepro_in_neg_table.size-1];
    } else {
    f -= i;
    f = tiltdrivepro_in_neg_table.data[i]*(1-f) + tiltdrivepro_in_neg_table.data[i+1]*f;
    }
    return copysign(f, x);
}


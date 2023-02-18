import("stdfaust.lib");

process = fi.conv(fcoeff);

// 128 tap 0~20000 -0.01dB // 30000~96000 -151.58dB
fcoeff = (0.200000003,
0.288000017,
0.332159996,
-0.449740797,
-0.283178478,
-0.16342932,
0.240353137,
-0.41975674,
0.155493334,
0.0965716988,
-0.174218506,
0.239860371,
-0.00281400839,
0.0510703847,
0.168975979,
0.0204230584,
0.130352408,
0.0969458893,
0.0651003644,
0.114502527,
0.0625997037,
0.0764238387,
0.0731616467,
0.0470526591,
0.0579451993,
0.0388748273,
0.0322872028,
0.0305418763,
0.0162300617,
0.0152165042,
0.00799270067,
0.00153163006,
-1.40097573e-05,
-0.00596009754,
-0.00787918642,
-0.00994974934,
-0.0126295742,
-0.012849981,
-0.0141066117,
-0.0143387746,
-0.0139447004,
-0.0139183998,
-0.0129893292,
-0.0121967439,
-0.0112597253,
-0.00999523699,
-0.00892739929,
-0.00766506512,
-0.0064438819,
-0.00531965215,
-0.00415151473,
-0.0031326164,
-0.00217067846,
-0.00128929946,
-0.000543868286,
0.000131288136,
0.000686885498,
0.00114130671,
0.00151168357,
0.00177911855,
0.00197251653,
0.00208992092,
0.0021380221,
0.00213410938,
0.00207799906,
0.00198350288,
0.00185823196,
0.00170723256,
0.00154107576,
0.00136366789,
0.00118142623,
0.000999963144,
0.000822166156,
0.000652579765,
0.000493461033,
0.000346784131,
0.000214474669,
9.70076653e-05,
-4.91265655e-06,
-9.12891555e-05,
-0.000162603988,
-0.000219318623,
-0.000262437825,
-0.000292988785,
-0.000312126387,
-0.000321192492,
-0.000321448373,
-0.000314238307,
-0.000300863467,
-0.000282541761,
-0.000260464818,
-0.000235698259,
-0.000209220612,
-0.000181909214,
-0.000154512731,
-0.000127682724,
-0.000101952442,
-7.77486639e-05,
-5.54033977e-05,
-3.51500203e-05,
-1.71418578e-05,
-1.45654781e-06,
1.18951857e-05,
2.29559319e-05,
3.18165221e-05,
3.86042884e-05,
4.3474909e-05,
4.66056881e-05,
4.818659e-05,
4.84151606e-05,
4.74903463e-05,
4.56074849e-05,
4.29548418e-05,
3.97098884e-05,
3.60369922e-05,
3.20857434e-05,
2.79893902e-05,
2.38645553e-05,
1.98108792e-05,
1.59112642e-05,
1.2232631e-05,
8.82660333e-06,
5.73075431e-06,
2.96982671e-06,
5.57043109e-07,
-1.50436949e-06,
-3.21981747e-06,
-4.60188767e-06,
-5.66895642e-06,
-6.44392048e-06,
-6.95296558e-06,
-7.22446703e-06,
-7.28801615e-06,
-7.17353032e-06,
-6.91053583e-06,
-6.52753079e-06,
-6.05147989e-06,
-5.50743698e-06,
-4.91823721e-06,
-4.30431282e-06,
-3.68359315e-06,
-3.07143819e-06,
-2.48071296e-06,
-1.92183097e-06,
-1.40290047e-06,
-9.29887221e-07,
-5.06788354e-07,
-1.35849433e-07,
1.822321e-07,
4.48083483e-07,
6.63445292e-07,
8.30962904e-07,
9.53985932e-07,
1.03638968e-06,
1.08240056e-06,
1.09644941e-06,
1.08303675e-06,
1.04661785e-06,
9.91507818e-07,
9.21801984e-07,
8.41314545e-07,
7.5353347e-07,
6.61588672e-07,
5.68233304e-07,
4.75838618e-07,
3.8639422e-07,
3.01522448e-07,
2.22494435e-07,
1.5025428e-07,
8.54471693e-08,
2.84485342e-08,
-2.06029647e-08,
-6.17733207e-08,
-9.53008907e-08,
-1.21564483e-07,
-1.41053704e-07,
-1.54340825e-07,
-1.62055002e-07,
-1.64859259e-07,
-1.63430016e-07,
-1.58439349e-07,
-1.5054033e-07,
-1.40354501e-07,
-1.28462517e-07,
-1.1539678e-07,
-1.01636452e-07,
-8.7604441e-08,
-7.36660439e-08,
-6.01292243e-08,
-4.72463135e-08,
-3.52162601e-08,
-2.41883917e-08,
-1.42662628e-08,
-5.51229951e-09,
2.04755612e-09,
8.41884873e-09,
1.36337253e-08,
1.7746272e-08,
2.08278284e-08,
2.29628156e-08,
2.42448017e-08,
2.47729055e-08,
2.46487524e-08,
2.39736657e-08,
2.2846411e-08,
2.13613216e-08,
1.96067198e-08,
1.76638526e-08,
1.56060445e-08,
1.34982043e-08,
1.13966019e-08,
9.34886746e-09,
7.39419415e-09,
5.56370816e-09,
3.8809711e-09,
2.36256836e-09,
1.01878006e-09,
-1.4568953e-10,
-1.13100262e-09,
-1.94143235e-09,
-2.58464672e-09,
-3.07102543e-09,
-3.41300321e-09,
-3.62447428e-09,
-3.72025299e-09,
-3.71559317e-09,
-3.62576169e-09,
-3.46569462e-09,
-3.24969274e-09,
-2.99119374e-09,
-2.7025866e-09,
-2.39508613e-09,
-2.07865236e-09,
-1.76194737e-09,
-1.45233048e-09,
-1.15588827e-09,
-8.7748353e-10,
-6.2082639e-10,
-3.88567095e-10,
-1.82389173e-10,
-3.12107284e-12,
1.49156312e-10,
2.7499869e-10,
3.75487641e-10,
4.52128612e-10,
5.06750197e-10,
5.4141408e-10,
5.58332214e-10,
5.59792102e-10,
5.48094126e-10,
5.25494537e-10,
4.94161601e-10,
4.56137933e-10,
4.13312051e-10,
3.67399111e-10,
3.19925947e-10,
2.72225381e-10,
2.25434268e-10,
1.80496604e-10,
1.38171252e-10,
9.90420732e-11,
6.35308611e-11,
3.19121916e-11,
4.32918371e-12,
-1.91898754e-11,
-3.87149479e-11,
-5.4397719e-11,
-6.64557506e-11,
-7.51575746e-11,
-8.08087486e-11,
-8.37391684e-11,
-8.42918443e-11,
-8.28129926e-11,
-7.96434169e-11,
-7.51116461e-11,
-6.95280222e-11,
-6.31804054e-11,
-5.63308566e-11,
-4.92134181e-11,
-4.20330125e-11,
-3.49651419e-11,
-2.81561441e-11,
-2.17243064e-11,
-1.5761354e-11,
-1.03342977e-11,
-5.48768583e-12,
-1.24589477e-12,
2.38434164e-12,
5.411443e-12,
7.85650191e-12,
9.75086852e-12,
1.11339037e-11,
1.20508273e-11,
1.25508128e-11,
1.26852556e-11,
1.25062183e-11,
1.20651726e-11,
1.14118628e-11,
1.05934246e-11,
9.65371619e-12,
8.63274556e-12,
7.56639022e-12,
6.48613932e-12,
5.4190619e-12,
4.38782864e-12,
3.41085883e-12,
2.50254743e-12,
1.67351636e-12,
9.30959303e-13,
2.78985954e-13,
-2.81014063e-13,
-7.49979614e-13,
-1.13081582e-12,
-1.42802588e-12,
-1.64736345e-12,
-1.79551686e-12,
-1.87980729e-12,
-1.90792716e-12,
-1.88770679e-12,
-1.82691058e-12,
-1.73306974e-12,
-1.61334054e-12,
-1.47439938e-12,
-1.3223565e-12,
-1.16270557e-12,
-1.00028243e-12,
-8.39261495e-13,
-6.83151777e-13,
-5.34819693e-13,
-3.96519704e-13,
-2.69933571e-13,
-1.56218205e-13,
-5.60593408e-14,
3.02750378e-14,
1.02876563e-13,
1.62140307e-13,
2.08710138e-13,
2.43426102e-13,
2.67275189e-13,
2.8134694e-13,
2.86791776e-13,
2.84786571e-13,
2.7650202e-13,
2.63077917e-13,
2.45599954e-13,
2.25084315e-13,
2.02463778e-13,
1.78578682e-13,
1.54172139e-13,
1.29886892e-13,
1.06265813e-13,
8.37547127e-14,
6.27061472e-14,
4.33859152e-14,
2.59795711e-14,
1.06004767e-14,
-2.7019307e-15,
-1.3933753e-14,
-2.31479603e-14,
-3.04361808e-14,
-3.59207293e-14,
-3.97471295e-14,
-4.20772629e-14,
-4.3083162e-14,
-4.29414838e-14,
-4.18287434e-14,
-3.99172543e-14,
-3.73717746e-14,
-3.43468235e-14,
-3.09847498e-14,
-2.74140655e-14,
-2.37488557e-14,
-2.00880877e-14,
-1.65156551e-14,
-1.31008443e-14,
-9.89875482e-15,
-6.95126651e-15,
-4.28810006e-15,
-1.92785292e-15,
1.20639117e-16,
1.85709581e-15,
3.28849869e-15,
4.42780671e-15,
5.29278956e-15,
5.90486912e-15,
6.28807444e-15,
6.46810367e-15,
6.4714664e-15,
6.32475351e-15,
6.05401552e-15,
5.6842255e-15,
5.23889285e-15,
4.73970584e-15,
4.20633428e-15,
3.65626769e-15,
3.10473701e-15,
2.56470925e-15,
2.04693262e-15,
1.56001075e-15,
1.11054151e-15,
7.03259057e-16,
3.41206418e-16,
2.59247404e-17,
-2.42360849e-16,
-4.64544103e-16,
-6.42449662e-16,
-7.78650177e-16,
-8.76293859e-16,
-9.38944281e-16,
-9.70436119e-16,
-9.74745822e-16,
-9.55876469e-16,
-9.1776401e-16,
-8.64193835e-16,
-7.98737829e-16,
-7.24702914e-16,
-6.45095528e-16,
-5.62596478e-16,
-4.79548974e-16,
-3.97955481e-16,
-3.19481875e-16,
-2.45471233e-16,
-1.76959913e-16,
-1.14701509e-16,
-5.91908013e-17,
-1.06919894e-17,
3.07323207e-17,
6.51929363e-17,
9.29446077e-17,
1.14358633e-16,
1.29897042e-16,
1.40087523e-16,
1.45501996e-16,
1.46736308e-16,
1.44392926e-16,
1.39066174e-16,
1.31329534e-16,
1.21725503e-16,
1.10758134e-16,
9.88865306e-17,
8.65216759e-17,
7.40238173e-17,
6.17017842e-17,
4.98139282e-17,
3.85693806e-17,
2.81309171e-17,
1.86181465e-17,
1.0111115e-17,
2.65463786e-18,
-3.73752544e-18,
-9.07836016e-18,
-1.34031053e-17,
-1.67652024e-17,
-1.92322967e-17,
-2.08824475e-17,
-2.18008645e-17,
-2.20767855e-17,
-2.18008149e-17,
-2.10626554e-17,
-1.99491363e-17,
-1.8542653e-17,
-1.69199497e-17,
-1.51511563e-17,
-1.32992196e-17,
-1.14194906e-17,
-9.5595973e-18,
-7.75954393e-18,
-6.05189366e-18,
-4.46216212e-18,
-3.00929237e-18,
-1.70620167e-18,
-5.60393358e-19,
4.2539164e-19,
1.25253064e-18,
1.92585369e-18,
2.45302189e-18,
2.84391312e-18,
3.11006737e-18,
3.26416451e-18,
3.31955728e-18,
3.2898677e-18,
3.18862575e-18,
3.02897328e-18,
2.82341621e-18,
2.58363826e-18,
2.3203446e-18,
2.04317127e-18,
1.7606194e-18,
1.48003146e-18,
1.20759529e-18,
9.48377151e-19,
7.06372797e-19,
4.84575689e-19,
2.85062828e-19,
1.09080178e-19,
-4.28567704e-20,
-1.70866421e-19,
-2.75602621e-19,
-3.58157754e-19,
-4.19972832e-19,
-4.62749753e-19,
-4.88374154e-19,
-4.98842952e-19,
-4.96202333e-19,
-4.82492332e-19,
-4.59700611e-19,
-4.29725369e-19,
-3.9434339e-19,
-3.55188744e-19,
-3.13736216e-19,
-2.71290167e-19,
-2.28982076e-19,
-1.87768655e-19,
-1.48437251e-19,
-1.11612824e-19,
-7.77678183e-20,
-4.72348467e-20,
-2.02195221e-20,
3.18489739e-21,
2.29824346e-20,
3.92601988e-20,
5.21732539e-20,
6.19310965e-20,
6.87845787e-20,
7.30136946e-20,
7.4916908e-20,
7.48013353e-20,
7.29744553e-20,
6.97367863e-20,
6.53761357e-20,
6.01628493e-20,
5.43460894e-20,
4.81514081e-20,
4.1779051e-20,
3.5403062e-20,
2.91713999e-20,
2.3206271e-20,
1.76053147e-20,
1.24430082e-20,
7.77242168e-21,
3.62729604e-21,
2.40827449e-23,
-3.03573447e-21,
-5.5634722e-21,
-7.58104084e-21,
-9.11881888e-21,
-1.02137031e-20,
-1.09072423e-20,
-1.12440034e-20,
-1.12700531e-20,
-1.10317031e-20,
-1.05743716e-20,
-9.9416883e-21,
-9.1747505e-21,
-8.31154068e-21,
-7.38655566e-21,
-6.43048916e-21,
-5.47013125e-21,
-4.52833039e-21,
-3.62404669e-21,
-2.77251595e-21,
-1.98545925e-21,
-1.27132874e-21,
-6.35623526e-22,
-8.11949539e-23,
3.91413066e-22,
7.83628119e-22,
1.09852312e-21,
1.34049243e-21,
1.51495201e-21,
1.62806219e-21,
1.6864692e-21,
1.6970822e-21,
1.66687051e-21,
1.60269535e-21,
1.51116719e-21,
1.39853088e-21,
1.27057557e-21,
1.13257182e-21,
9.89225719e-22,
8.44656318e-22,
7.02390239e-22,
5.65367649e-22,
4.35963462e-22,
3.16018486e-22,
2.06877144e-22,
1.09430558e-22,
2.41655816e-23,
-4.87860204e-23,
-1.09597354e-22,
-1.58695363e-22,
-1.96713788e-22,
-2.24446503e-22,
-2.42805202e-22,
-2.52780948e-22,
-2.55408422e-22,
-2.5173637e-22,
-2.42800633e-22,
-2.2960222e-22,
-2.13089885e-22,
-1.94145371e-22,
-1.73573671e-22,
-1.52096018e-22,
-1.30345406e-22,
-1.08865872e-22,
-8.81128286e-23,
-6.84562513e-23,
-5.01849737e-23,
-3.35120545e-23,
-1.85816186e-23,
-5.47566143e-24,
5.77836117e-24,
1.5199895e-23,
2.28478447e-23,
2.88131951e-23,
3.32120791e-23,
3.61793068e-23,
3.78624394e-23,
3.84164091e-23,
3.79988808e-23,
3.67661783e-23,
3.48698293e-23,
3.24538575e-23,
2.96525525e-23,
2.65887935e-23,
2.33730572e-23,
2.01025772e-23,
1.68612045e-23,
1.37194681e-23,
1.07349375e-23,
7.95284915e-24,
5.40695125e-24,
3.12042254e-24,
1.10695309e-24,
-6.28141934e-25,
-2.08678539e-24,
-3.27700117e-24,
-4.21180213e-24,
-4.90814109e-24,
-5.38592495e-24,
-5.66710298e-24,
-5.77485823e-24,
-5.73288018e-24,
-5.56474079e-24,
-5.2933754e-24,
-4.94063506e-24,
-4.52695837e-24,
-4.07110839e-24,
-3.58998213e-24,
-3.09851626e-24,
-2.60963411e-24,
-2.13423853e-24,
-1.68128485e-24,
-1.25785311e-24,
-8.69271272e-25,
-5.19260838e-25,
-2.10087538e-25,
5.72678874e-26,
2.82939212e-25,
4.68003466e-25,
6.14313893e-25,
7.24339035e-25,
8.01013948e-25,
8.47600869e-25,
8.67565755e-25,
8.64465532e-25,
8.41854116e-25,
8.03199734e-25,
7.51817032e-25,
6.90815503e-25,
6.23057134e-25,
5.51127712e-25,
4.7731947e-25,
4.03621545e-25,
3.31721188e-25,
2.63008198e-25,
1.98589656e-25,
1.39306328e-25,
8.57526995e-26,
3.83019025e-26,
-2.87091867e-27,
-3.77618717e-26,
-6.65128257e-26,
-8.93863037e-26,
-1.0674088e-25,
-1.19008555e-25,
-1.26673867e-25,
-1.30254321e-25,
-1.30283842e-25,
-1.27297868e-25,
-1.21820573e-25,
-1.14354645e-25,
-1.05372897e-25,
-9.53116537e-26,
-8.45665115e-26,
-7.34889631e-26,
-6.23852774e-26,
-5.15159698e-26,
-4.10968983e-26,
-3.13009037e-26,
-2.22603035e-26,
-1.40700153e-26,
-6.79097073e-27,
-4.5383711e-28,
4.93712873e-27,
9.40014898e-27,
1.29721736e-26,
1.57051513e-26,
1.76625756e-26,
1.89162944e-26,
1.95435636e-26,
1.96244481e-26,
1.92395764e-26,
1.84681428e-26,
1.7386331e-26,
1.60660136e-26,
1.45736937e-26,
1.29698378e-26,
1.13083504e-26,
9.63632037e-27,
7.99400133e-27,
6.41485125e-27,
4.92583161e-27,
3.5477594e-27,
2.2957239e-27,
1.17964162e-27,
2.04782927e-28,
-6.27644295e-28,
-1.31989892e-27,
-1.87714518e-27,
-2.30688447e-27,
-2.61843272e-27,
-2.82243339e-27,
-2.93040622e-27,
-2.95435304e-27,
-2.90640432e-27,
-2.79852219e-27,
-2.642248e-27,
-2.44849714e-27,
-2.22740808e-27,
-1.98821452e-27,
-1.73918003e-27,
-1.48754486e-27,
-1.23951638e-27,
-1.00028449e-27,
-7.74048867e-28,
-5.64076521e-28,
-3.72765715e-28,
-2.01719615e-28,
-5.18335216e-29,
7.66237608e-29,
1.83918088e-28,
2.70763939e-28,
3.38242072e-28,
3.87715963e-28,
4.20760371e-28,
4.39092094e-28,
4.44507953e-28,
4.38833056e-28,
4.23872721e-28,
4.01374249e-28,
3.72995888e-28,
3.40280958e-28,
3.04640228e-28,
2.6733925e-28,
2.2949049e-28,
1.92051688e-28,
1.55825969e-28,
1.21467583e-28,
8.9488708e-29,
6.02691489e-29,
3.40677658e-29,
1.10343771e-29,
-8.77683741e-30,
-2.53944666e-29,
-3.89165803e-29,
-4.94979102e-29,
-5.73379016e-29,
-6.26691766e-29,
-6.57472279e-29,
-6.68410899e-29,
-6.6225015e-29,
-6.41713478e-29,
-6.09446086e-29,
-5.67964565e-29,
-5.19618427e-29,
-4.66561585e-29,
-4.10731132e-29,
-3.53836303e-29,
-2.9735281e-29,
-2.42523646e-29,
-1.90366544e-29,
-1.41683624e-29,
-9.70752124e-30,
-5.69574823e-30,
-2.15796952e-30,
8.9562833e-31,
3.46753718e-30,
5.57104953e-30,
7.22825883e-30,
8.4682267e-30,
9.32530013e-30,
9.83745982e-30,
1.00449434e-29,
9.98895001e-30,
9.71054706e-30,
9.24974424e-30,
8.64473366e-30,
7.93125185e-30,
7.14218255e-30,
6.30716365e-30,
5.45242691e-30,
4.60071411e-30,
3.77124242e-30,
2.97983248e-30,
2.23902653e-30,
1.55830358e-30,
9.44332352e-31,
4.01219873e-31,
-6.9175421e-32,
-4.66959441e-31,
-7.93902459e-31,
-1.05313987e-30,
-1.24890266e-30,
-1.38624883e-30,
-1.47082188e-30,
-1.50863528e-30,
-1.50587259e-30,
-1.46872339e-30,
-1.40323839e-30,
-1.3152075e-30,
-1.21006874e-30,
-1.09283669e-30,
-9.68043578e-31,
-8.39716552e-31,
-7.1135539e-31,
-5.85930472e-31,
-4.65898581e-31,
-3.53219055e-31,
-2.49386275e-31,
-1.55464759e-31,
-7.21286096e-32,
2.93915212e-34,
6.1776601e-32,
1.12549987e-31,
1.53057287e-31,
1.83912133e-31,
2.05858495e-31,
2.19734125e-31,
2.26436582e-31,
2.26893896e-31,
2.22038282e-31,
2.12783991e-31,
2.00008953e-31,
1.84539941e-31,
1.67141073e-31,
1.48505914e-31,
1.29251705e-31,
1.09916891e-31,
9.09606345e-32,
7.27636939e-32,
5.56321568e-32,
3.98011014e-32,
2.543999e-32,
1.26589765e-32,
1.51479326e-33,
-7.98204378e-33,
-1.58607249e-32,
-2.21834798e-32,
-2.70390179e-32,
-3.05366543e-32,
-3.28005212e-32,
-3.39645486e-32,
-3.41678827e-32,
-3.35508363e-32,
-3.2251536e-32,
-3.04029741e-32,
-2.81308258e-32,
-2.55515267e-32,
-2.27710697e-32,
-1.98840747e-32,
-1.69733479e-32,
-1.41097555e-32,
-1.13523462e-32,
-8.74883858e-33,
-6.33615994e-33,
-4.14125448e-33,
-2.18199853e-33,
-4.68077396e-34,
9.97924633e-34,
2.21954658e-33,
3.20545922e-33,
3.96844946e-33,
4.52453984e-33,
4.89211288e-33,
5.09113069e-33,
5.14243918e-33,
5.06715244e-33,
4.88612153e-33,
4.61949003e-33,
4.28634721e-33,
3.90443277e-33,
3.48993216e-33,
3.0573495e-33,
2.61940806e-33,
2.1870414e-33,
1.76939983e-33,
1.37391063e-33,
1.0063735e-33,
6.71061302e-34,
3.70857953e-34,
1.07404356e-34,
-1.18761124e-34,
-3.08039397e-34,
-4.61623887e-34,
-5.81354828e-34,
-6.69574347e-34,
-7.29001051e-34,
-7.6260848e-34,
-7.73516976e-34,
-7.64901521e-34,
-7.39908537e-34,
-7.01587834e-34,
-6.52837374e-34,
-5.96358131e-34,
-5.34622315e-34,
-4.69850143e-34,
-4.0399725e-34,
-3.38748337e-34,
-2.75520631e-34,
-2.15470089e-34,
-1.59505033e-34,
-1.08302129e-34,
-6.23258469e-35,
-2.18498184e-35,
1.30208772e-35,
4.23264894e-35,
6.62298313e-35,
8.49939985e-35,
9.89611271e-35,
1.08532578e-34,
1.14150764e-34,
1.16282576e-34,
1.15405099e-34,
1.1199286e-34,
1.06507411e-34,
9.93884721e-35,
9.10471975e-35,
8.18608871e-35,
7.21694011e-35,
6.22730186e-35,
5.24313643e-35,
4.28636496e-35,
3.37496596e-35,
2.5231507e-35,
1.74161312e-35,
1.03780754e-35,
4.15883466e-36,
-1.2096576e-36,
-5.74537895e-36,
-9.46420033e-36,
-1.23977331e-35,
-1.46073511e-35,
-1.61430573e-35,
-1.70734183e-35,
-1.74709197e-35,
-1.74024104e-35,
-1.69432444e-35,
-1.6161776e-35,
-1.51240636e-35,
-1.38942844e-35,
-1.25284969e-35,
-1.10794669e-35,
-9.59330023e-36,
-8.10947107e-36,
-6.66246579e-36,
-5.27981854e-36,
-3.98383997e-36,
-2.79153703e-36,
-1.71460405e-36,
-7.60679154e-37,
6.68246335e-38,
7.67894429e-37,
1.34534069e-36,
1.8045501e-36,
2.1527225e-36,
2.39857204e-36,
2.55188782e-36,
2.62306804e-36,
2.62288383e-36,
2.5621099e-36,
2.45129082e-36,
2.30055415e-36,
2.11939734e-36,
1.91660842e-36,
1.70013781e-36,
1.47329947e-36,
1.25350547e-36,
1.03202632e-36,
8.35885072e-37,
6.26748069e-37,
4.52820832e-37,
2.90700712e-37,
1.37651544e-37,
1.75613498e-38,
-9.31470218e-38,
-1.83368603e-37,
-2.55508749e-37,
-3.12914745e-37,
-3.50669403e-37,
-3.74898212e-37,
-3.91221097e-37,
-3.90172791e-37,
-3.8383955e-37,
-3.69601864e-37,
-3.46572454e-37,
-3.21729652e-37,
-2.91612363e-37,
-2.59435994e-37,
-2.20869906e-37,
-1.93076463e-37,
-1.56172651e-37,
-1.26157679e-37,
-9.70814419e-38,
-7.20257097e-38,
-4.18243815e-38,
-2.7094784e-38,
0,
2.06254585e-38,
3.61299614e-38,
5.56228834e-38,
4.54192978e-38,
5.28723699e-38,
6.73720704e-38,
6.26510342e-38,
7.34737891e-38,
6.06567567e-38,
8.38406961e-38,
6.24955573e-38,
6.51334456e-38,
6.9009241e-38,
4.01827295e-38,
5.40878002e-38,
4.81495681e-38,
3.78293945e-38,
1.92069619e-38,
3.09810891e-38,
2.28133156e-38,
0,
1.30052885e-38,
0);

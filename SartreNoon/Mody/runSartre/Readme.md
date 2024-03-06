# runSartre
**runSartre** folder contains the files for merging different methods of **runB** Folder.

### test7.C
This file is used for extracting vector meson rapidity data from **example_Pb.root** file.

### test8.C
This file is used for extracting the energy of thevirtual photon exchanged between the two UPC ions.

### test10.C
This file is used for caluclating the cross-section *d sigma/dy* from the rapidity of the vector meson.

## extract.C
This File contains the merging of **test7**, **test8** and **test10**.

Hence this file can extract the vector meson *rapidity*, *photon flux* and calculate *cross-section* at once.

The output is **extract.root**.

The output contains 3 histograms 

1) The 1st histogram is **Rapidity**
2) The 2nd histogram is **PhotonK**
3) The 3rd histogram is **XSection**

## run.C
The file **run.C** contains the merged code of **extract.C** and **runB.C** from *runB* folder.

It can extract data from **example_Pb.root** and calculate the breakup probabilities.

It has 3 main methods:
1) **extract()**
2) **runSartre()**
3) **computerModelBreakup()**

## runSartre()
This is still in the building stage.

In the final satage I expect it to run Sartre and compute the neutron generation all in one place.

As of now it generated events and check neutron breakup probabilities and store everything in the **noon.root** file. which is it's output file.

It takes input as the files **Rapidity** and **PhotonK**.

## computeModelBreakup()
This file taked the *cross-section* as input and guve us neutron generation/ nuclear breakup probabilities inthe form of 3 graphs. 

It takes input as **XSection** i.e the *cross-section* "*d sigma/dy*" and the *photon flux* **PhotonK**.

From this it calculates and plots nuclear breakups.

### My thoughts

1) I think we don't need 2 methods **runSartre()** and **computeModelbreaup()**. In the final stage bith methods can be merged into a single one, at least thaw's what I am assuming. 

2) Need more understanidn of the two methods.
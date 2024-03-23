# Info
1) I have created the program *run.C* which will create the histogram of the rapidity values from the *example_Pb.root* File
2) As of now I am storing it in a file *temp.root*. But it will be later on deleted or will be moved into *example_Pb.root*.
3) The method of working is as follows

> I will create the hitogram of rapidity
> Take random points from the hostgram
> Genrate the phootn flux corresponding to that rapidity value 
> Then I will generate the events i.e neutron generation
> It will then store it into *noon1.root* for now

4) The furthur plan is to takeevent by event. For that I will do the following

> I have set the number random picking of the rapidity to be sameas that of the number of events generated in sartre.
> I will Generate a event corresponding to a event in sartre.
> Then I will write the generated neutron data to the corresponding event to sartre.
> In this will randomly events will be allocated to each events of sartre.

{
	struct BeamHit
	{
	    int id;
	    double q;
    	    pair<double, double> hit;
	};
	
	
	vector<pair<int, int>> hitDat;
	
	hitDat = {
		make_pair(1, 1),
		make_pair(1, 5),
		make_pair(2, 1), 
		make_pair(2, 5), 
		make_pair(3, 2), 
		make_pair(3, 4), 
		make_pair(4, 3), 
		make_pair(9, 9), 
		make_pair(9, 8),
		make_pair(8, 7) 
	};
	
	vector<vector<BeamHit>> branch; // [cluster][hit]
	//找cluster
	for (int i = 0; i < (int)hitDat.size(); i++)
	{
	        double xhit = hitDat[i].first;
	        double yhit = hitDat[i].second;
	
		BeamHit _hit_; 
        	_hit_.hit = hitDat[i];

		int checked = false;

	        for (int ii = 0; ii < (int)branch.size(); ii++)
	        {
		    int slen = branch[ii].size();
	            for (int jj = 0; jj < slen; jj++)
		    {
	                if (fabs(xhit - branch[ii][jj].hit.first) <= 1 && fabs(yhit - branch[ii][jj].hit.second) <= 1)
	                {
	                    branch[ii].push_back(_hit_);
			    checked = true;
	                }
	            }
                }

		if(!checked)
		{
			branch.resize(branch.size()+1);
			branch[branch.size()-1].push_back(_hit_);
		}
	}
	
	cout<<branch.size()<<endl;
	for(int ii=0; ii<(int)branch.size(); ii++)
 	{
		cout<<"\n"<<ii<<":"<<endl;
		for(int jj=0; jj<(int)branch[ii].size(); jj++)
		{
			cout<<"("<<branch[ii][jj].hit.first<<", "<<branch[ii][jj].hit.second<<") ";
		}
		cout<<endl;
	}


	if(branch.size()<=1) return;

    for (int i = 0; i < branch.size() - 1; i++)
        for (int j = 1; j < branch.size(); j++)
        {
            int mergeFlag = 0;

            for (int ii = 0; ii < (int)branch[i].size(); ii++)
                for (int jj = 0; jj < (int)branch[j].size(); jj++)
                    if (fabs(branch[i][ii].hit.first - branch[j][jj].hit.first) < 1 &&
                        fabs(branch[i][ii].hit.second - branch[j][jj].hit.second) < 1)
                    {
                        mergeFlag = 1;
                        break;
                    }

            if (mergeFlag) //cluster_i和_j需要合并
            {
                for (int jj = 0; jj < (int)branch[j].size(); jj++)
                {
                    int duplicateFlag = 0;
                    for (int ii = 0; ii < (int)branch[i].size(); ii++)
                        if (branch[i][ii].hit.first == branch[j][jj].hit.first &&
                            branch[i][ii].hit.second == branch[j][jj].hit.second)
                        {
                            duplicateFlag = 1;
                            break;
                        }

                    if(duplicateFlag == 0)
                        branch[i].push_back(branch[j][jj]);
                }

                branch[j].resize(0);
            }
        }


	for(int ii=(int)branch.size(); ii>=0; ii--)
	   if(branch[ii].size()==0) branch.erase(branch.begin()+ii);

        cout<<"\nafter mearge: "<<branch.size()<<endl;
        for(int ii=0; ii<(int)branch.size(); ii++)
        {
                cout<<ii<<":"<<endl;
                for(int jj=0; jj<(int)branch[ii].size(); jj++)
                {
                        cout<<"("<<branch[ii][jj].hit.first<<", "<<branch[ii][jj].hit.second<<") ";
                }
                cout<<endl;
        }



}


{
    //LiF
    //ref[eng_] := Sqrt[ 1 + (34.76 * (12.632^2 - eng^2))/((12.632^2 - eng^2)^2 + 0.33^2*eng^2) + 236.6/(18.37^2 - eng^2)]
    {
        for (int wave = 120; wave < 221; wave++)
        {
            double eng = 1239.8521 / wave;
            double eng2 = eng * eng;
            double tmp = 12.632 * 12.632 - eng2;
            double ref = 1 + 34.76 * tmp / (tmp * tmp + 0.33 * 0.33 * eng2) + 236.6 / (18.37 * 18.37 - eng2);
            ref = sqrt(ref);
            cout << wave << ", " << ref << ", " << endl;
        }
    }
}
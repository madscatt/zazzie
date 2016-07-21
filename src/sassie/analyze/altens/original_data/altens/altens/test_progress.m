function test_progress(k,steps);

percentage=(k/steps)*100;

if ((percentage==10)|(percentage==20)|(percentage==30)|(percentage==40)|(percentage==50)|...
        (percentage==60)|(percentage==70)|(percentage==80)|(percentage==90)|(percentage==100)),
    
    fprintf(' %2.0f%%',percentage);
end
return

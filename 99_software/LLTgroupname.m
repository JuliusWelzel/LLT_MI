function g_name = LLTgroupname(ID)

SUBJ = str2double(extractAfter({ID},'SUBJ'));

    if SUBJ <70 & SUBJ >= 30
        g_name = "Old";

    elseif SUBJ >= 70;
        g_name = "Young";

    elseif SUBJ < 30;
        g_name = "Stroke";
    
    end
end


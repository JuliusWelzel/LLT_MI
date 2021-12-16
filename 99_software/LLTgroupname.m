function g_name = LLTgroupname(ID)

if contains(ID,'SUBJ')
    SUBJ = str2double(extractAfter({ID},'SUBJ'));
        if SUBJ <70 & SUBJ >= 30
            g_name = "old_nic";
        elseif SUBJ >= 70;
            g_name = "young_nic";
        elseif SUBJ < 30;
            g_name = "stroke_nic";
        end
else
    g_name = "old_juw";
end
end


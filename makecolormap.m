% ***************************************************************
% *** Source Code is mainly written for research purposes. The codes are
% *** having copyrights and required proper citations whenever it is used.
% *** Originated by:
% ***       Mr. Arka Roy (email: arka.phy@gmail.com)
% ***       Solid Earth Research Group, National Centre for Earth Science Studies,
% ***       Ministry of Earth Sciences, Government of India
% ***       Thiruvanthapuram, Kerala, India
% ****************************************************************

%%Matlab function for creating colormap as per user choice 
%% all color name and corresponding color representation can be found from 
% the web address  http://kodama.fubuki.info/wiki/wiki.cgi/GrADS/script/color.gs?lang=en

function ccmap=makecolormap(colors, values)
    sz1=size(colors);
    sz2=values+1;
    N=sz2;
    N1=sz2;
    c(1:sz1(1,2)-1)=0;
    while N>0
        for k=1:sz1(1,2)-1
            if N>0
                c(k)=c(k)+1;
            else
                c(k)=c(k)+0;
            end
            N=N1-sum(c(:));
        end
    end
    for i=1:sz1(1,2)
     cols=colors(i);
     t=strcmp(cols,'black');if (t==1) ; r=0 ; g=0 ; b=0 ; end
     t=strcmp(cols,'navy');if (t==1) ; r=0 ; g=0 ; b=128 ; end
     t=strcmp(cols,'darkblue');if (t==1) ; r=0 ; g=0 ; b=139 ; end
     t=strcmp(cols,'mediumblue');if (t==1) ; r=0 ; g=0 ; b=205 ; end
     t=strcmp(cols,'blue');if (t==1) ; r=0 ; g=0 ; b=255 ; end
     t=strcmp(cols,'darkgreen');if (t==1) ; r=0 ; g=100 ; b=0 ; end
     t=strcmp(cols,'green');if (t==1) ; r=0 ; g=128 ; b=0 ; end
     t=strcmp(cols,'teal');if (t==1) ; r=0 ; g=128 ; b=128 ; end
     t=strcmp(cols,'darkcyan');if (t==1) ; r=0 ; g=139 ; b=139 ; end
     t=strcmp(cols,'deepskyblue');if (t==1) ; r=0 ; g=191 ; b=255 ; end
     t=strcmp(cols,'darkturquoise');if (t==1) ; r=0 ; g=206 ; b=209 ; end
     t=strcmp(cols,'mediumspringgreen');if (t==1) ; r=0 ; g=250 ; b=154 ; end
     t=strcmp(cols,'lime');if (t==1) ; r=0 ; g=255 ; b=0 ; end
     t=strcmp(cols,'springgreen');if (t==1) ; r=0 ; g=255 ; b=127 ; end
     t=strcmp(cols,'aqua');if (t==1) ; r=0 ; g=255 ; b=255 ; end
     t=strcmp(cols,'cyan');if (t==1) ; r=0 ; g=255 ; b=255 ; end
     t=strcmp(cols,'midnightblue');if (t==1) ; r=25 ; g=25 ; b=112 ; end
     t=strcmp(cols,'dodgerblue');if (t==1) ; r=30 ; g=144 ; b=255 ; end
     t=strcmp(cols,'lightseagreen');if (t==1) ; r=32 ; g=178 ; b=170 ; end
     t=strcmp(cols,'forestgreen');if (t==1) ; r=34 ; g=139 ; b=34 ; end
     t=strcmp(cols,'seagreen');if (t==1) ; r=46 ; g=139 ; b=87 ; end
     t=strcmp(cols,'darkslategray');if (t==1) ; r=47 ; g=79 ; b=79 ; end
     t=strcmp(cols,'limegreen');if (t==1) ; r=50 ; g=205 ; b=50 ; end
     t=strcmp(cols,'mediumseagreen');if (t==1) ; r=60 ; g=179 ; b=113 ; end
     t=strcmp(cols,'turquoise');if (t==1) ; r=64 ; g=224 ; b=208 ; end
     t=strcmp(cols,'royalblue');if (t==1) ; r=65 ; g=105 ; b=225 ; end
     t=strcmp(cols,'steelblue');if (t==1) ; r=70 ; g=130 ; b=180 ; end
     t=strcmp(cols,'darkslateblue');if (t==1) ; r=72 ; g=61 ; b=139 ; end
     t=strcmp(cols,'mediumturquoise');if (t==1) ; r=72 ; g=209 ; b=204 ; end
     t=strcmp(cols,'indigo');if (t==1) ; r=75 ; g=0 ; b=130 ; end
     t=strcmp(cols,'darkolivegreen');if (t==1) ; r=85 ; g=107 ; b=47 ; end
     t=strcmp(cols,'cadetblue');if (t==1) ; r=95 ; g=158 ; b=160 ; end
     t=strcmp(cols,'cornflowerblue');if (t==1) ; r=100 ; g=149 ; b=237 ; end
     t=strcmp(cols,'mediumaquamarine');if (t==1) ; r=102 ; g=205 ; b=170 ; end
     t=strcmp(cols,'dimgray');if (t==1) ; r=105 ; g=105 ; b=105 ; end
     t=strcmp(cols,'slateblue');if (t==1) ; r=106 ; g=90 ; b=205 ; end
     t=strcmp(cols,'olivedrab');if (t==1) ; r=107 ; g=142 ; b=35 ; end
     t=strcmp(cols,'slategray');if (t==1) ; r=112 ; g=128 ; b=144 ; end
     t=strcmp(cols,'lightslategray');if (t==1) ; r=119 ; g=136 ; b=153 ; end
     t=strcmp(cols,'mediumslateblue');if (t==1) ; r=123 ; g=104 ; b=238 ; end
     t=strcmp(cols,'lawngreen');if (t==1) ; r=124 ; g=252 ; b=0 ; end
     t=strcmp(cols,'chartreuse');if (t==1) ; r=127 ; g=255 ; b=0 ; end
     t=strcmp(cols,'aquamarine');if (t==1) ; r=127 ; g=255 ; b=212 ; end
     t=strcmp(cols,'maroon');if (t==1) ; r=128 ; g=0 ; b=0 ; end
     t=strcmp(cols,'purple');if (t==1) ; r=128 ; g=0 ; b=128 ; end
     t=strcmp(cols,'olive');if (t==1) ; r=128 ; g=128 ; b=0 ; end
     t=strcmp(cols,'gray');if (t==1) ; r=128 ; g=128 ; b=128 ; end
     t=strcmp(cols,'skyblue');if (t==1) ; r=135 ; g=206 ; b=235 ; end
     t=strcmp(cols,'lightskyblue');if (t==1) ; r=135 ; g=206 ; b=250 ; end
     t=strcmp(cols,'blueviolet');if (t==1) ; r=138 ; g=43 ; b=226 ; end
     t=strcmp(cols,'darkred');if (t==1) ; r=139 ; g=0 ; b=0 ; end
     t=strcmp(cols,'darkmagenta');if (t==1) ; r=139 ; g=0 ; b=139 ; end
     t=strcmp(cols,'saddlebrown');if (t==1) ; r=139 ; g=69 ; b=19 ; end
     t=strcmp(cols,'darkseagreen');if (t==1) ; r=143 ; g=188 ; b=143 ; end
     t=strcmp(cols,'lightgreen');if (t==1) ; r=144 ; g=238 ; b=144 ; end
     t=strcmp(cols,'mediumpurple');if (t==1) ; r=147 ; g=112 ; b=219 ; end
     t=strcmp(cols,'darkviolet');if (t==1) ; r=148 ; g=0 ; b=211 ; end
     t=strcmp(cols,'palegreen');if (t==1) ; r=152 ; g=251 ; b=152 ; end
     t=strcmp(cols,'darkorchid');if (t==1) ; r=153 ; g=50 ; b=204 ; end
     t=strcmp(cols,'yellowgreen');if (t==1) ; r=154 ; g=205 ; b=50 ; end
     t=strcmp(cols,'sienna');if (t==1) ; r=160 ; g=82 ; b=45 ; end
     t=strcmp(cols,'brown');if (t==1) ; r=165 ; g=42 ; b=42 ; end
     t=strcmp(cols,'darkgray');if (t==1) ; r=169 ; g=169 ; b=169 ; end
     t=strcmp(cols,'lightblue');if (t==1) ; r=173 ; g=216 ; b=230 ; end
     t=strcmp(cols,'greenyellow');if (t==1) ; r=173 ; g=255 ; b=47 ; end
     t=strcmp(cols,'paleturquoise');if (t==1) ; r=175 ; g=238 ; b=238 ; end
     t=strcmp(cols,'lightsteelblue');if (t==1) ; r=176 ; g=196 ; b=222 ; end
     t=strcmp(cols,'powderblue');if (t==1) ; r=176 ; g=224 ; b=230 ; end
     t=strcmp(cols,'firebrick');if (t==1) ; r=178 ; g=34 ; b=34 ; end
     t=strcmp(cols,'darkgoldenrod');if (t==1) ; r=184 ; g=134 ; b=11 ; end
     t=strcmp(cols,'mediumorchid');if (t==1) ; r=186 ; g=85 ; b=211 ; end
     t=strcmp(cols,'rosybrown');if (t==1) ; r=188 ; g=143 ; b=143 ; end
     t=strcmp(cols,'darkkhaki');if (t==1) ; r=189 ; g=183 ; b=107 ; end
     t=strcmp(cols,'silver');if (t==1) ; r=192 ; g=192 ; b=192 ; end
     t=strcmp(cols,'mediumvioletred');if (t==1) ; r=199 ; g=21 ; b=133 ; end
     t=strcmp(cols,'indianred');if (t==1) ; r=205 ; g=92 ; b=92 ; end
     t=strcmp(cols,'peru');if (t==1) ; r=205 ; g=133 ; b=63 ; end
     t=strcmp(cols,'chocolate');if (t==1) ; r=210 ; g=105 ; b=30 ; end
     t=strcmp(cols,'tan');if (t==1) ; r=210 ; g=180 ; b=140 ; end
     t=strcmp(cols,'lightgray');if (t==1) ; r=211 ; g=211 ; b=211 ; end
     t=strcmp(cols,'thistle');if (t==1) ; r=216 ; g=191 ; b=216 ; end
     t=strcmp(cols,'orchid');if (t==1) ; r=218 ; g=112 ; b=214 ; end
     t=strcmp(cols,'goldenrod');if (t==1) ; r=218 ; g=165 ; b=32 ; end
     t=strcmp(cols,'palevioletred');if (t==1) ; r=219 ; g=112 ; b=147 ; end
     t=strcmp(cols,'crimson');if (t==1) ; r=220 ; g=20 ; b=60 ; end
     t=strcmp(cols,'gainsboro');if (t==1) ; r=220 ; g=220 ; b=220 ; end
     t=strcmp(cols,'plum');if (t==1) ; r=221 ; g=160 ; b=221 ; end
     t=strcmp(cols,'burlywood');if (t==1) ; r=222 ; g=184 ; b=135 ; end
     t=strcmp(cols,'lightcyan');if (t==1) ; r=224 ; g=255 ; b=255 ; end
     t=strcmp(cols,'lavender');if (t==1) ; r=230 ; g=230 ; b=250 ; end
     t=strcmp(cols,'darksalmon');if (t==1) ; r=233 ; g=150 ; b=122 ; end
     t=strcmp(cols,'violet');if (t==1) ; r=238 ; g=130 ; b=238 ; end
     t=strcmp(cols,'palegoldenrod');if (t==1) ; r=238 ; g=232 ; b=170 ; end
     t=strcmp(cols,'lightcoral');if (t==1) ; r=240 ; g=128 ; b=128 ; end
     t=strcmp(cols,'khaki');if (t==1) ; r=240 ; g=230 ; b=140 ; end
     t=strcmp(cols,'aliceblue');if (t==1) ; r=240 ; g=248 ; b=255 ; end
     t=strcmp(cols,'honeydew');if (t==1) ; r=240 ; g=255 ; b=240 ; end
     t=strcmp(cols,'azure');if (t==1) ; r=240 ; g=255 ; b=255 ; end
     t=strcmp(cols,'sandybrown');if (t==1) ; r=244 ; g=164 ; b=96 ; end
     t=strcmp(cols,'wheat');if (t==1) ; r=245 ; g=222 ; b=179 ; end
     t=strcmp(cols,'beige');if (t==1) ; r=245 ; g=245 ; b=220 ; end
     t=strcmp(cols,'whitesmoke');if (t==1) ; r=245 ; g=245 ; b=245 ; end
     t=strcmp(cols,'mintcream');if (t==1) ; r=245 ; g=255 ; b=250 ; end
     t=strcmp(cols,'ghostwhite');if (t==1) ; r=248 ; g=248 ; b=255 ; end
     t=strcmp(cols,'salmon');if (t==1) ; r=250 ; g=128 ; b=114 ; end
     t=strcmp(cols,'antiquewhite');if (t==1) ; r=250 ; g=235 ; b=215 ; end
     t=strcmp(cols,'linen');if (t==1) ; r=250 ; g=240 ; b=230 ; end
     t=strcmp(cols,'lightgoldenrodyellow');if (t==1) ; r=250 ; g=250 ; b=210 ; end
     t=strcmp(cols,'oldlace');if (t==1) ; r=253 ; g=245 ; b=230 ; end
     t=strcmp(cols,'red');if (t==1) ; r=255 ; g=0 ; b=0 ; end
     t=strcmp(cols,'fuchsia');if (t==1) ; r=255 ; g=0 ; b=255 ; end
     t=strcmp(cols,'magenta');if (t==1) ; r=255 ; g=0 ; b=255 ; end
     t=strcmp(cols,'deeppink');if (t==1) ; r=255 ; g=20 ; b=147 ; end
     t=strcmp(cols,'orangered');if (t==1) ; r=255 ; g=69 ; b=0 ; end
     t=strcmp(cols,'tomato');if (t==1) ; r=255 ; g=99 ; b=71 ; end
     t=strcmp(cols,'hotpink');if (t==1) ; r=255 ; g=105 ; b=180 ; end
     t=strcmp(cols,'coral');if (t==1) ; r=255 ; g=127 ; b=80 ; end
     t=strcmp(cols,'darkorange');if (t==1) ; r=255 ; g=140 ; b=0 ; end
     t=strcmp(cols,'lightsalmon');if (t==1) ; r=255 ; g=160 ; b=122 ; end
     t=strcmp(cols,'orange');if (t==1) ; r=255 ; g=165 ; b=0 ; end
     t=strcmp(cols,'lightpink');if (t==1) ; r=255 ; g=182 ; b=193 ; end
     t=strcmp(cols,'pink');if (t==1) ; r=255 ; g=192 ; b=203 ; end
     t=strcmp(cols,'gold');if (t==1) ; r=255 ; g=215 ; b=0 ; end
     t=strcmp(cols,'peachpuff');if (t==1) ; r=255 ; g=218 ; b=185 ; end
     t=strcmp(cols,'navajowhite');if (t==1) ; r=255 ; g=222 ; b=173 ; end
     t=strcmp(cols,'moccasin');if (t==1) ; r=255 ; g=228 ; b=181 ; end
     t=strcmp(cols,'bisque');if (t==1) ; r=255 ; g=228 ; b=196 ; end
     t=strcmp(cols,'mistyrose');if (t==1) ; r=255 ; g=228 ; b=225 ; end
     t=strcmp(cols,'blanchedalmond');if (t==1) ; r=255 ; g=235 ; b=205 ; end
     t=strcmp(cols,'papayawhip');if (t==1) ; r=255 ; g=239 ; b=213 ; end
     t=strcmp(cols,'lavenderblush');if (t==1) ; r=255 ; g=240 ; b=245 ; end
     t=strcmp(cols,'seashell');if (t==1) ; r=255 ; g=245 ; b=238 ; end
     t=strcmp(cols,'cornsilk');if (t==1) ; r=255 ; g=248 ; b=220 ; end
     t=strcmp(cols,'lemonchiffon');if (t==1) ; r=255 ; g=250 ; b=205 ; end
     t=strcmp(cols,'floralwhite');if (t==1) ; r=255 ; g=250 ; b=240 ; end
     t=strcmp(cols,'snow');if (t==1) ; r=255 ; g=250 ; b=250 ; end
     t=strcmp(cols,'yellow');if (t==1) ; r=255 ; g=255 ; b=0 ; end
     t=strcmp(cols,'lightyellow');if (t==1) ; r=255 ; g=255 ; b=224 ; end
     t=strcmp(cols,'ivory');if (t==1) ; r=255 ; g=255 ; b=240 ; end
     t=strcmp(cols,'white');if (t==1) ; r=255 ; g=255 ; b=255 ; end
     ccolors(i,1:3)=[r/255 g/255 b/255];
    end
    z=1;
    mm=0;
    for i=1:sz1(1,2)-1

        mm=mm+c(i);
        if i==1
            cmap=interpmap(ccolors(i,1:3),ccolors(i+1,1:3),c(i));
            ccmap(z:mm,1:3)=cmap;
        else
            cmap=interpmap(ccolors(i,1:3),ccolors(i+1,1:3),c(i)+1);
            ccmap(z:mm,1:3)=cmap(2:end,1:3);
        end
        z=z+c(i);

    end

    function cmap=interpmap(colorstart, colorend, n)
            for ii=1:3
                cmap(1:n,ii)=linspace(colorstart(ii),colorend(ii),n);
            end
    end
end
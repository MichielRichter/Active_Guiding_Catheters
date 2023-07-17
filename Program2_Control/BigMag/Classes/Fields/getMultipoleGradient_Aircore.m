function grad_B = getMultipoleGradient_Aircore(in1)
%GETMULTIPOLEGRADIENT_AIRCORE
%    GRAD_B = GETMULTIPOLEGRADIENT_AIRCORE(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    17-Jun-2020 19:16:19

x = in1(1,:);
y = in1(2,:);
z = in1(3,:);
t2 = conj(x);
t3 = conj(y);
t4 = conj(z);
t5 = x.^2;
t6 = y.^2;
t7 = z.^2;
t8 = t2.^2;
t9 = t2.^3;
t11 = t2.^5;
t12 = t3.^2;
t13 = t4.^2;
t14 = t5+t6+t7;
t10 = t8.^2;
t15 = t14.^(5.0./2.0);
t16 = t14.^(7.0./2.0);
t17 = t14.^(9.0./2.0);
t18 = t14.^(1.1e+1./2.0);
t19 = t14.^(1.3e+1./2.0);
t20 = conj(t15);
t21 = conj(t16);
t22 = conj(t17);
t23 = conj(t18);
t24 = conj(t19);
t25 = 1.0./t20;
t26 = 1.0./t21;
t27 = 1.0./t22;
t28 = 1.0./t23;
t29 = 1.0./t24;
t30 = t2.*t25.*1.926794026857413e-6;
t31 = t3.*t25.*1.926794026857413e-6;
t32 = t4.*t25.*1.926794026857413e-6;
t33 = t25.*5.612494655481929e-8;
t34 = t2.*t3.*t4.*t26.*9.633970134287063e-6;
t35 = t3.*t8.*t26.*9.633970134287063e-6;
t36 = t4.*t8.*t26.*9.633970134287063e-6;
t40 = t3.*t4.*t26.*2.806247327740964e-7;
t41 = t2.*t3.*t26.*8.418741983222893e-7;
t42 = t2.*t4.*t26.*8.418741983222893e-7;
t43 = t2.*t26.*7.85740607956154e-9;
t44 = t3.*t26.*7.85740607956154e-9;
t45 = t4.*t26.*7.85740607956154e-9;
t46 = t8.*t26.*2.806247327740964e-7;
t47 = t26.*2.332698908387632e-11;
t56 = t9.*t27.*1.833394751897693e-8;
t57 = t3.*t9.*t27.*1.964373129418675e-6;
t58 = t4.*t9.*t27.*1.964373129418675e-6;
t59 = t2.*t3.*t4.*t27.*5.500184255693078e-8;
t60 = t3.*t10.*t28.*1.650055276707923e-7;
t61 = t4.*t10.*t28.*1.650055276707923e-7;
t62 = t8.*t27.*3.265778471742684e-10;
t63 = t3.*t4.*t8.*t27.*1.964373129418675e-6;
t64 = t3.*t8.*t27.*1.100036851138616e-7;
t65 = t4.*t8.*t27.*1.100036851138616e-7;
t66 = t3.*t4.*t27.*1.632889235871342e-10;
t69 = t10.*t28.*4.898667707614027e-10;
t70 = t2.*t3.*t27.*8.164446179356711e-10;
t71 = t2.*t4.*t27.*8.164446179356711e-10;
t72 = t3.*t4.*t9.*t28.*1.650055276707923e-7;
t75 = t3.*t9.*t28.*4.898667707614027e-9;
t76 = t4.*t9.*t28.*4.898667707614027e-9;
t77 = t3.*t11.*t29.*5.388534478375429e-9;
t78 = t4.*t11.*t29.*5.388534478375429e-9;
t81 = t3.*t4.*t8.*t28.*2.939200624568416e-9;
t82 = t3.*t4.*t10.*t29.*5.388534478375429e-9;
t37 = -t34;
t38 = -t35;
t39 = -t36;
t48 = -t40;
t49 = -t41;
t50 = -t42;
t51 = -t43;
t52 = -t44;
t53 = -t45;
t54 = -t46;
t55 = -t47;
t67 = -t60;
t68 = -t61;
t73 = -t69;
t74 = -t72;
t79 = -t75;
t80 = -t76;
t83 = -t81;
t84 = t37+t48+t59+t63+t66+t74+t82+t83;
t85 = t31+t38+t49+t52+t57+t64+t67+t70+t77+t79;
t86 = t32+t39+t50+t53+t58+t65+t68+t71+t78+t80;
grad_B = reshape([t25.*1.683748396644579e-7-t26.*1.166349454193816e-10+t2.*t25.*5.780382080572238e-6-t2.*t26.*3.92870303978077e-8-t8.*t26.*1.683748396644579e-6+t8.*t27.*2.449333853807013e-9-t9.*t26.*9.633970134287063e-6+t9.*t27.*1.833394751897693e-7+t10.*t27.*1.964373129418675e-6-t10.*t28.*7.34800156142104e-9-t11.*t28.*1.650055276707923e-7+t8.^3.*t29.*5.388534478375429e-9,t85,t86,t85,t30+t33+t51+t54+t55+t56+t62+t73-t12.*t26.*2.806247327740964e-7+t12.*t27.*1.632889235871342e-10-t2.*t12.*t26.*9.633970134287063e-6+t2.*t12.*t27.*5.500184255693078e-8+t8.*t12.*t27.*1.964373129418675e-6-t8.*t12.*t28.*2.939200624568416e-9-t9.*t12.*t28.*1.650055276707923e-7+t10.*t12.*t29.*5.388534478375429e-9,t84,t86,t84,t30+t33+t51+t54+t55+t56+t62+t73-t13.*t26.*2.806247327740964e-7+t13.*t27.*1.632889235871342e-10-t2.*t13.*t26.*9.633970134287063e-6+t2.*t13.*t27.*5.500184255693078e-8+t8.*t13.*t27.*1.964373129418675e-6-t8.*t13.*t28.*2.939200624568416e-9-t9.*t13.*t28.*1.650055276707923e-7+t10.*t13.*t29.*5.388534478375429e-9],[3,3]);
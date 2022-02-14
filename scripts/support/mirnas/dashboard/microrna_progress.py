"""
Store the updated and newly committed microRNA families.
"""

updated_families = []
new_commits = []

# Rfam 14.3
updated_families += ['RF00725', 'RF00821', 'RF00680', 'RF00823', 'RF00709', 'RF00740', 'RF00986', 'RF00765', 'RF00748', 'RF00666', 'RF00760', 'RF00967', 'RF00786', 'RF00681', 'RF00971', 'RF00970', 'RF01899', 'RF00820', 'RF00799', 'RF01901', 'RF00790', 'RF00736', 'RF00363', 'RF00787', 'RF00739', 'RF00813', 'RF00845', 'RF00718', 'RF00968', 'RF00777', 'RF00670', 'RF00758', 'RF00847', 'RF00868', 'RF00966', 'RF00723', 'RF00864', 'RF00867', 'RF00866']
new_commits = ['RF03555', 'RF03554', 'RF03553', 'RF03552', 'RF03550', 'RF03549', 'RF03548', 'RF03527', 'RF03526', 'RF03525', 'RF03524', 'RF03523', 'RF03522', 'RF03521', 'RF03520', 'RF03519', 'RF03518', 'RF03517', 'RF03516', 'RF03514', 'RF03512', 'RF03511', 'RF03510', 'RF03509', 'RF03508', 'RF03507', 'RF03506', 'RF03505', 'RF03504', 'RF03503', 'RF03502', 'RF03501', 'RF03500', 'RF03499', 'RF03498', 'RF03497', 'RF03496', 'RF03495', 'RF03494', 'RF03493', 'RF03492', 'RF03491', 'RF03490', 'RF03489', 'RF03488', 'RF03487', 'RF03486', 'RF03485', 'RF03484', 'RF03483', 'RF03482', 'RF03481', 'RF03480', 'RF03479', 'RF03478', 'RF03477', 'RF03476', 'RF03475', 'RF03474', 'RF03473', 'RF03472', 'RF03471', 'RF03470', 'RF03469', 'RF03468', 'RF03467', 'RF03467', 'RF03466', 'RF03465', 'RF03464', 'RF03463', 'RF03462', 'RF03461', 'RF03460', 'RF03458', 'RF03457', 'RF03456', 'RF03455', 'RF03454', 'RF03453', 'RF03452', 'RF03451', 'RF03450', 'RF03449', 'RF03448', 'RF03447', 'RF03446', 'RF03445', 'RF03444', 'RF03443', 'RF03442', 'RF03441', 'RF03440', 'RF03439', 'RF03438', 'RF03437', 'RF03436', 'RF03435', 'RF03434', 'RF03433', 'RF03432', 'RF03431', 'RF03430', 'RF03429', 'RF03428', 'RF03427', 'RF03426', 'RF03425', 'RF03424', 'RF03423', 'RF03422', 'RF03421', 'RF03420', 'RF03419', 'RF03418', 'RF03417', 'RF03416', 'RF03414', 'RF03413', 'RF03412', 'RF03411', 'RF03410', 'RF03409', 'RF03408', 'RF03407', 'RF03406', 'RF03405', 'RF03404', 'RF03403', 'RF03402', 'RF03401', 'RF03400', 'RF03399', 'RF03398', 'RF03397', 'RF03396', 'RF03395', 'RF03394', 'RF03393', 'RF03392', 'RF03391', 'RF03390', 'RF03389', 'RF03388', 'RF03387', 'RF03386', 'RF03385', 'RF03384', 'RF03383', 'RF03382', 'RF03381', 'RF03380', 'RF03379', 'RF03378', 'RF03377', 'RF03376', 'RF03375', 'RF03374', 'RF03373', 'RF03372', 'RF03371', 'RF03370', 'RF03369', 'RF03368', 'RF03367', 'RF03366', 'RF03365', 'RF03364', 'RF03363', 'RF03362', 'RF03361', 'RF03360', 'RF03359', 'RF03358', 'RF03357', 'RF03356', 'RF03355', 'RF03354', 'RF03353', 'RF03352', 'RF03351', 'RF03350', 'RF03349', 'RF03348', 'RF03347', 'RF03346', 'RF03345', 'RF03344', 'RF03343', 'RF03342', 'RF03341', 'RF03340', 'RF03339', 'RF03338', 'RF03337', 'RF03336', 'RF03335', 'RF03334', 'RF03333', 'RF03332', 'RF03331', 'RF03330', 'RF03329', 'RF03328', 'RF03327', 'RF03326', 'RF03325', 'RF03324', 'RF03323', 'RF03322', 'RF03321', 'RF03320', 'RF03319', 'RF03318', 'RF03317', 'RF03316', 'RF03315', 'RF03314', 'RF03313', 'RF03312', 'RF03311', 'RF03310', 'RF03309', 'RF03308', 'RF03307', 'RF03306', 'RF03305', 'RF03304', 'RF03303', 'RF03302', 'RF03301', 'RF03300', 'RF03299', 'RF03298', 'RF03297', 'RF03296', 'RF03295', 'RF03294', 'RF03293', 'RF03292', 'RF03291', 'RF03290', 'RF03289', 'RF03288', 'RF03287', 'RF03286', 'RF03285', 'RF03284', 'RF03283', 'RF03281', 'RF03280', 'RF03279', 'RF03278', 'RF03277', 'RF03276', 'RF03275', 'RF03274', 'RF03272', 'RF03271', 'RF03270', 'RF03269', 'RF03268', 'RF03267', 'RF03266', 'RF03265', 'RF03264', 'RF03263', 'RF03262', 'RF03261', 'RF03260', 'RF03259', 'RF03258', 'RF03257', 'RF03256', 'RF03255', 'RF03254', 'RF03253', 'RF03252', 'RF03251', 'RF03249', 'RF03248', 'RF03247', 'RF03246', 'RF03244', 'RF03243', 'RF03242', 'RF03241', 'RF03240', 'RF03239', 'RF03238', 'RF03237', 'RF03236', 'RF03235', 'RF03234', 'RF03233', 'RF03232', 'RF03231', 'RF03230', 'RF03229', 'RF03228', 'RF03227', 'RF03226', 'RF03225', 'RF03224', 'RF03223', 'RF03222', 'RF03221', 'RF03220', 'RF03219', 'RF03218', 'RF03217', 'RF03216', 'RF03215', 'RF03214', 'RF03213', 'RF03212', 'RF03211', 'RF03210', 'RF03209', 'RF03208', 'RF03207', 'RF03206', 'RF03205', 'RF03204', 'RF03203', 'RF03202', 'RF03201', 'RF03200', 'RF03199', 'RF03198', 'RF03197', 'RF03196', 'RF03195', 'RF03194', 'RF03193', 'RF03192', 'RF03191', 'RF03190', 'RF03189', 'RF03188', 'RF03187', 'RF03186', 'RF03185', 'RF03184', 'RF03183', 'RF03182', 'RF03181', 'RF03180', 'RF03179', 'RF03178', 'RF03177', 'RF03176', 'RF03175', 'RF03174', 'RF03173', 'RF03172', 'RF03171']

# Rfam 14.4
new_commits += ['RF04051', 'RF04050', 'RF04049', 'RF04048', 'RF04047', 'RF04046', 'RF04045', 'RF04044', 'RF04043', 'RF04042', 'RF04041', 'RF04040', 'RF04039', 'RF04038', 'RF04037', 'RF04036', 'RF04035', 'RF04034', 'RF04033', 'RF04032', 'RF04031', 'RF04030', 'RF04029', 'RF04028', 'RF04027', 'RF04026', 'RF04025', 'RF04024', 'RF04023', 'RF04022', 'RF04021', 'RF04020', 'RF04019', 'RF04018', 'RF04017', 'RF04016', 'RF04015', 'RF04014', 'RF04013', 'RF04012', 'RF04011', 'RF04010', 'RF04009', 'RF04008', 'RF04007', 'RF04006', 'RF04005', 'RF04004', 'RF04003', 'RF04002', 'RF04001', 'RF04000', 'RF03999', 'RF03998', 'RF03997', 'RF03996', 'RF03995', 'RF03994', 'RF03993', 'RF03992', 'RF03991', 'RF03990', 'RF03989', 'RF03988', 'RF03987', 'RF03986', 'RF03985', 'RF03984', 'RF03983', 'RF03982', 'RF03981', 'RF03980', 'RF03979', 'RF03978', 'RF03977', 'RF03976', 'RF03975', 'RF03974', 'RF03973', 'RF03972', 'RF03971', 'RF03970', 'RF03969', 'RF03968', 'RF03967', 'RF03966', 'RF03965', 'RF03964', 'RF03963', 'RF03962', 'RF03961', 'RF03960', 'RF03959', 'RF03958', 'RF03957', 'RF03956', 'RF03955', 'RF03954', 'RF03953', 'RF03952', 'RF03951', 'RF03950', 'RF03949', 'RF03948', 'RF03947', 'RF03946', 'RF03945', 'RF03944', 'RF03943', 'RF03942', 'RF03941', 'RF03940', 'RF03939', 'RF03938', 'RF03937', 'RF03936', 'RF03935', 'RF03934', 'RF03933', 'RF03932', 'RF03931', 'RF03930', 'RF03929', 'RF03928', 'RF03927', 'RF03926', 'RF03925', 'RF03924', 'RF03923', 'RF03922', 'RF03921', 'RF03920', 'RF03919', 'RF03918', 'RF03917', 'RF03916', 'RF03915', 'RF03914', 'RF03913', 'RF03912', 'RF03911', 'RF03910', 'RF03909', 'RF03908', 'RF03907', 'RF03906', 'RF03905', 'RF03904', 'RF03903', 'RF03902', 'RF03901', 'RF03900', 'RF03899', 'RF03898', 'RF03897', 'RF03896', 'RF03895', 'RF03894', 'RF03893', 'RF03892', 'RF03891', 'RF03890', 'RF03889', 'RF03888', 'RF03887', 'RF03886', 'RF03885', 'RF03884', 'RF03883', 'RF03882', 'RF03881', 'RF03880', 'RF03879', 'RF03878', 'RF03877', 'RF03876', 'RF03875', 'RF03874', 'RF03873', 'RF03872', 'RF03871', 'RF03870', 'RF03869', 'RF03868', 'RF03867', 'RF03866', 'RF03865', 'RF03864', 'RF03863', 'RF03862', 'RF03861', 'RF03860', 'RF03859', 'RF03858', 'RF03857', 'RF03856', 'RF03855', 'RF03854', 'RF03853', 'RF03852', 'RF03851', 'RF03850', 'RF03849', 'RF03848', 'RF03847', 'RF03846', 'RF03845', 'RF03844', 'RF03843', 'RF03842', 'RF03841', 'RF03840', 'RF03839', 'RF03838', 'RF03837', 'RF03836', 'RF03835', 'RF03834', 'RF03833', 'RF03832', 'RF03831', 'RF03830', 'RF03829', 'RF03828', 'RF03827', 'RF03826', 'RF03825', 'RF03824', 'RF03823', 'RF03822', 'RF03821', 'RF03820', 'RF03819', 'RF03818', 'RF03817', 'RF03816', 'RF03815', 'RF03814', 'RF03813', 'RF03812', 'RF03811', 'RF03810', 'RF03809', 'RF03808', 'RF03807', 'RF03806', 'RF03805', 'RF03804', 'RF03803', 'RF03802', 'RF03801', 'RF03800', 'RF03799', 'RF03798', 'RF03797', 'RF03796', 'RF03795', 'RF03794', 'RF03793', 'RF03792', 'RF03791', 'RF03790', 'RF03789', 'RF03788', 'RF03787', 'RF03786', 'RF03785', 'RF03784', 'RF03783', 'RF03782', 'RF03781', 'RF03780', 'RF03779', 'RF03778', 'RF03777', 'RF03776', 'RF03775', 'RF03774', 'RF03773', 'RF03772', 'RF03771', 'RF03770', 'RF03769', 'RF03768', 'RF03767', 'RF03766', 'RF03765', 'RF03764', 'RF03763', 'RF03762', 'RF03761', 'RF03760', 'RF03759', 'RF03758', 'RF03757', 'RF03756', 'RF03755', 'RF03754', 'RF03753', 'RF03752', 'RF03751', 'RF03750', 'RF03749', 'RF03748', 'RF03747', 'RF03746', 'RF03745', 'RF03744', 'RF03743', 'RF03742', 'RF03741', 'RF03740', 'RF03739', 'RF03738', 'RF03737', 'RF03736', 'RF03735', 'RF03734', 'RF03733', 'RF03732', 'RF03731', 'RF03730', 'RF03729', 'RF03728', 'RF03727', 'RF03726', 'RF03725', 'RF03724', 'RF03723', 'RF03722', 'RF03721', 'RF03720', 'RF03719', 'RF03718', 'RF03717', 'RF03716', 'RF03715', 'RF03714', 'RF03713', 'RF03712', 'RF03711', 'RF03710', 'RF03709', 'RF03708', 'RF03707', 'RF03706', 'RF03705', 'RF03704', 'RF03703', 'RF03702', 'RF03701', 'RF03700', 'RF03699', 'RF03698', 'RF03697', 'RF03696', 'RF03695', 'RF03694', 'RF03693', 'RF03692', 'RF03691', 'RF03690', 'RF03689', 'RF03688', 'RF03687', 'RF03686', 'RF03685', 'RF03684', 'RF03683', 'RF03682', 'RF03681', 'RF03680', 'RF03679', 'RF03678', 'RF03677', 'RF03676', 'RF03675', 'RF03674', 'RF03673', 'RF03672', 'RF03671', 'RF03670', 'RF03669', 'RF03668', 'RF03667', 'RF03666', 'RF03665', 'RF03664', 'RF03663', 'RF03662', 'RF03661', 'RF03660', 'RF03659', 'RF03658', 'RF03657', 'RF03656', 'RF03655', 'RF03654', 'RF03653', 'RF03652', 'RF03651', 'RF03650', 'RF03649', 'RF03648', 'RF03647', 'RF03646', 'RF03645', 'RF03644', 'RF03643', 'RF03642', 'RF03641', 'RF03640', 'RF03639', 'RF03638', 'RF03637', 'RF03636', 'RF03635', 'RF03634', 'RF03633', 'RF03632', 'RF03631', 'RF03630', 'RF03629', 'RF03628', 'RF03627', 'RF03626', 'RF03625', 'RF03624', 'RF03623', 'RF03622', 'RF03621', 'RF03620', 'RF03619', 'RF03618', 'RF03617', 'RF03616', 'RF03615', 'RF03614', 'RF03613', 'RF03612', 'RF03611', 'RF03610', 'RF03609', 'RF03608', 'RF03607', 'RF03606', 'RF03605', 'RF03604', 'RF03603', 'RF03602', 'RF03601', 'RF03600', 'RF03599', 'RF03598', 'RF03597', 'RF03596', 'RF03595', 'RF03594', 'RF03593', 'RF03592', 'RF03591', 'RF03590', 'RF03589', 'RF03588', 'RF03587', 'RF03586', 'RF03585', 'RF03584', 'RF03583', 'RF03582', 'RF03581', 'RF03580', 'RF03579', 'RF03578', 'RF03577', 'RF03576', 'RF03575', 'RF03574', 'RF03573', 'RF03572', 'RF03571', 'RF03570', 'RF03569', 'RF03568', 'RF03567', 'RF03566', 'RF03565', 'RF03564', 'RF03563', 'RF03562', 'RF03561', 'RF03560', 'RF03559', 'RF03558', 'RF03557', 'RF03556']

# Rfam 14.5
updated_families += ['RF00253', 'RF00365', 'RF00366', 'RF00658', 'RF00664', 'RF00676', 'RF00682', 'RF00683', 'RF00687', 'RF00699', 'RF00707', 'RF00712', 'RF00715', 'RF00719', 'RF00730', 'RF00731', 'RF00733', 'RF00744', 'RF00745', 'RF00746', 'RF00752', 'RF00753', 'RF00755', 'RF00757', 'RF00762', 'RF00763', 'RF00767', 'RF00770', 'RF00772', 'RF00781', 'RF00783', 'RF00794', 'RF00798', 'RF00801', 'RF00807', 'RF00808', 'RF00811', 'RF00817', 'RF00822', 'RF00826', 'RF00836', 'RF00839', 'RF00849', 'RF00852', 'RF00853', 'RF00854', 'RF00857', 'RF00858', 'RF00859', 'RF00861', 'RF00862', 'RF00869', 'RF00870', 'RF00877', 'RF00892', 'RF00897', 'RF00900', 'RF00902', 'RF00904', 'RF00910', 'RF00912', 'RF00914', 'RF00919', 'RF00921', 'RF00922', 'RF00927', 'RF00931', 'RF00932', 'RF00936', 'RF00937', 'RF00939', 'RF00941', 'RF00947', 'RF00961', 'RF00962', 'RF00964', 'RF00979', 'RF00981', 'RF00985', 'RF00987', 'RF00988', 'RF00993', 'RF01000', 'RF01001', 'RF01010', 'RF01011', 'RF01012', 'RF01015', 'RF01035', 'RF01038', 'RF01044', 'RF01064', 'RF01916', 'RF01918', 'RF01920', 'RF01921', 'RF01938', 'RF01939', 'RF01943', 'RF01945', 'RF01997', 'RF02009', 'RF02016', 'RF02017', 'RF02018', 'RF02019', 'RF02023', 'RF02024', 'RF02026', 'RF02092', 'RF02097', 'RF02244']

# Rfam 14.6
new_commits += ['RF04052', 'RF04053', 'RF04054', 'RF04055', 'RF04056', 'RF04057', 'RF04059', 'RF04060', 'RF04061', 'RF04062', 'RF04063', 'RF04064', 'RF04065', 'RF04066', 'RF04067', 'RF04068', 'RF04069', 'RF04070', 'RF04071', 'RF04072', 'RF04073', 'RF04074', 'RF04075', 'RF04076', 'RF04077', 'RF04078', 'RF04079', 'RF04080', 'RF04081', 'RF04082', 'RF04083', 'RF04084', 'RF04085', 'RF04087', 'RF04088', 'RF04089', 'RF04090', 'RF04091', 'RF04092', 'RF04093', 'RF04094', 'RF04095', 'RF04096', 'RF04097', 'RF04098', 'RF04099', 'RF04100', 'RF04101', 'RF04102', 'RF04103', 'RF04104', 'RF04105', 'RF04106', 'RF04107', 'RF04108', 'RF04109', 'RF04110', 'RF04111', 'RF04112', 'RF04113', 'RF04114', 'RF04115', 'RF04116', 'RF04117', 'RF04118', 'RF04119', 'RF04120', 'RF04121', 'RF04122', 'RF04123', 'RF04124', 'RF04125', 'RF04126', 'RF04127', 'RF04128', 'RF04129', 'RF04130', 'RF04131', 'RF04132', 'RF04133', 'RF04134', 'RF04135', 'RF04136', 'RF04137', 'RF04138', 'RF04139', 'RF04140', 'RF04141', 'RF04142', 'RF04143', 'RF04144', 'RF04145', 'RF04146', 'RF04147', 'RF04148', 'RF04149', 'RF04150', 'RF04151', 'RF04152', 'RF04153', 'RF04154', 'RF04155', 'RF04156', 'RF04157', 'RF04158', 'RF04159', 'RF04160', 'RF04161', 'RF04162', 'RF04163', 'RF04164', 'RF04165', 'RF04166', 'RF04167', 'RF04168', 'RF04169', 'RF04170', 'RF04171', 'RF04172', 'RF04173', 'RF04174', 'RF04175', 'RF04176', 'RF04185', 'RF04186', 'RF04187']

# Rfam 14.7
updated_families += ['RF00129', 'RF00130', 'RF00143', 'RF00144', 'RF00178', 'RF00244', 'RF00249', 'RF00255', 'RF00637', 'RF00644', 'RF00646', 'RF00649', 'RF00650', 'RF00651', 'RF00653', 'RF00656', 'RF00659', 'RF00662', 'RF00663', 'RF00668', 'RF00669', 'RF00671', 'RF00673', 'RF00674', 'RF00675', 'RF00684', 'RF00686', 'RF00694', 'RF00696', 'RF00698', 'RF00701', 'RF00703', 'RF00705', 'RF00708', 'RF00710', 'RF00716', 'RF00720', 'RF00721', 'RF00727', 'RF00729', 'RF00732', 'RF00734', 'RF00735', 'RF00737', 'RF00749', 'RF00750', 'RF00751', 'RF00756', 'RF00764', 'RF00766', 'RF00769', 'RF00773', 'RF00774', 'RF00784', 'RF00785', 'RF00788', 'RF00792', 'RF00793', 'RF00795', 'RF00796', 'RF00803', 'RF00806', 'RF00809', 'RF00810', 'RF00812', 'RF00814', 'RF00816', 'RF00818', 'RF00819', 'RF00829', 'RF00832', 'RF00835', 'RF00843', 'RF00850', 'RF00856', 'RF00873', 'RF00878', 'RF00879', 'RF00891', 'RF00894', 'RF00901', 'RF00905', 'RF00935', 'RF00949', 'RF00950', 'RF00955', 'RF00963', 'RF00969', 'RF00974', 'RF00977', 'RF00984', 'RF00992', 'RF00996', 'RF00997', 'RF01007', 'RF01009', 'RF01013', 'RF01014', 'RF01024', 'RF01025', 'RF01026', 'RF01027', 'RF01030', 'RF01034', 'RF01040', 'RF01041', 'RF01314', 'RF01895', 'RF01896', 'RF01900', 'RF01912', 'RF01914', 'RF01915', 'RF01919', 'RF01923', 'RF01926', 'RF01941', 'RF01944', 'RF02008', 'RF02520', 'RF03459']

# Rfam 14.8
updated_families += []
new_commits += []

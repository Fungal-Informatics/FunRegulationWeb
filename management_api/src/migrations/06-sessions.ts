import { MigrationDefinitionWithName } from "umzug/lib/migrationsList";
import { zap } from "../helpers/sql";

const migration: MigrationDefinitionWithName = {
	name: "06-sessions",
	async up() {
		await zap.transaction.single/*sql*/`            
            create table 
                sessions
            ( 
                key text primary key,
                "userId" uuid references users(id) not null
            );
        `;
	},
	async down() {
		await zap.transaction.single/*sql*/`
            drop table sessions
        `;
	},
};

export default migration;
